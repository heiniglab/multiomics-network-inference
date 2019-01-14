#' -----------------------------------------------------------------------------
#' Perform LD based pruning for the eQTLGen trans-eQTL results
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Jan 11 23:18:07 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params and load data.")
# ------------------------------------------------------------------------------

# input
feqtl <- snakemake@input$eqtl
eqtl <- fread(feqtl)
eqtl_by_chromosome <- split(eqtl, f=eqtl$SNPChr)

fdosage <- snakemake@input$dosage
findividuals <- snakemake@input$individuals
individuals <- read.table(findividuals, header=F)[,1]

# output
fpruned <- snakemake@output[[1]]

# params
ld_max_dist <- as.numeric(snakemake@params$ld_max_dist)
r2_min <- as.numeric(snakemake@params$r2_min)

# ------------------------------------------------------------------------------
print("Pruning SNPs for all chromosomes.")
# ------------------------------------------------------------------------------
pruned <- lapply(names(eqtl_by_chromosome), function(chr) {
  print(paste0("Current chromosome: ", chr))

  # get eQTL subset for current chromosome
  eqtl_subs <- eqtl_by_chromosome[[chr]]
  eqtl_subs <- eqtl_subs[order(eqtl_subs$SNPPos, decreasing = F),]

  # get unique set of SNPs
  snps <- with(eqtl_subs, GRanges(paste0("chr", chr),
                                         IRanges(SNPPos, width=1)))
  names(snps) <- eqtl_subs$SNP
  snps <- unique(snps)

  # get KORA genotypes
  geno <- scan_snps(snps, fdosage, individuals)

  # get a list of MAFs to be used later on
  MAF <- rowMeans(geno) / 2
  MAF[MAF>0.5] <- 1 - MAF[MAF>0.5]

  # TODO cope with missing genotypes for some SNPs

  # calculate pairwise r^2
  pw_r2 <- cor(t(geno))^2

  # check each SNP pair (within distance threshold) for large r^2 and define
  # LD clusters

  # list of cluster and cluster identifier (id of first SNP)
  cluster <- list()
  clustern <- names(snps)[1]
  cluster[[clustern]] <- clustern

  for(i in 1:(length(snps)-1)) {
    # first SNP
    s1 <- snps[i]
    for (j in (i + 1):length(snps)) {
      # second SNP
      s2 <- snps[j]
      r2 <-pw_r2[names(s1), names(s2)]

      # check LD criteria distance and r^2
      d12 <- distance(s1,s2)
      if (r2 > 0.2 & d12 <= ld_max_dist) {
        # add to current cluster
        cluster[[clustern]] <- c(cluster[[clustern]], names(s2))
        # we can skip this SNP in the next round
        i <- j + 1
        # we replace the current s1 SNP (to also adjust distance)
        s1 <- s2

      } else {
        # in that case, we do not have to test the other SNPs since the list
        # is sorted by genomic position
        # we create a new cluster with the next SNP (i+1 is here always in range)
        clustern <- names(snps)[i+1]
        cluster[[clustern]] <- clustern
        break
      }
    }
  }

  # now we select a representative SNP for each cluster using
  # 1) the number of trans associations
  # 2) the MAF
  sentinels <- unlist(lapply(cluster, function(cl_snps) {
    # get the number of trans associations per SNP
    eqtl <- eqtl[eqtl$SNP %in% cl_snps,]
    counts <- tapply(eqtl$GeneSymbol, eqtl$SNP, length)
    counts <- counts[counts == max(counts)]
    # check whether we have to check MAF as well
    if(length(counts) > 1) {
      mafs <- MAF[names(counts)]
      # select one (random) snp with the highest MAF
      names(mafs[mafs == max(mafs)])[1]
    } else {
      names(counts)
    }
  }))

  # annotate the eqtl with the respective sentinel
  eqtl_subs <- cbind(eqtl_subs, sentinel=NA_character_)

  for(i in 1:length(cluster)) {
    # get all snps
    cl_snps <- cluster[[i]]

    # get the sentinel
    sentinel <- sentinels[names(cluster)[i]]

    # set the sentinel in the eQTL list
    eqtl_subs[eqtl_subs$SNP %in% cl_snps]$sentinel <- sentinel
  }
  eqtl_subs
})

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fpruned, pruned)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
