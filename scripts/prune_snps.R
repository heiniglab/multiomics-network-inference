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
library(parallel)
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params and load data.")
# ------------------------------------------------------------------------------

# params
ld_max_dist <- as.numeric(snakemake@params$ld_max_dist)
r2_min <- as.numeric(snakemake@params$r2_min)
threads <- snakemake@threads
fdr_cutoff <- as.numeric(snakemake@params$fdr_cutoff)
tempdir <- snakemake@params$tempdir
print(paste0("using ", tempdir, " as temp directory."))

# input
feqtl <- snakemake@input$eqtl
eqtl <- fread(paste0("zcat ", feqtl))
# we also filter to fdr < fdrcutoff
eqtl <- eqtl[eqtl$FDR < fdr_cutoff]

# we noticed that there are not only trans eQTLs in this table, so filter for
# distinct chromosome sbetween snp and gene
eqtl <- eqtl[eqtl$SNPChr != eqtl$GeneChr]

eqtl_by_chromosome <- split(eqtl, f=eqtl$SNPChr)

fdosage <- snakemake@input$dosage
findividuals <- snakemake@input$individuals
individuals <- read.table(findividuals, header=F)[,1]

# output
ffull <- snakemake@output$full
fpruned <- snakemake@output$pruned

# ------------------------------------------------------------------------------
print("Pruning SNPs for all chromosomes.")
# ------------------------------------------------------------------------------
full <- mclapply(names(eqtl_by_chromosome), function(chr) {
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
  geno <- scan_snps(snps, fdosage, individuals, tempdir=tempdir)

  # remove SNPs for which we do not have any genotypes
  snps <- snps[names(snps)[names(snps) %in% rownames(geno)]]
  eqtl_subs <- eqtl_subs[eqtl_subs$SNP %in% names(snps)]

  # get MAFs to be used later on
  MAF <- rowMeans(geno) / 2
  MAF[MAF>0.5] <- 1 - MAF[MAF>0.5]

  # calculate pairwise r^2
  pw_r2 <- cor(t(geno))^2

  # check each SNP pair (within distance threshold) for large r^2 and define
  # LD clusters

  # list of cluster and first cluster identifier (id of first SNP)
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
      if (r2 > r2_min & d12 <= ld_max_dist) {
        # add to current cluster
        cluster[[clustern]] <- c(cluster[[clustern]], names(s2))

        # we replace the current reference s1 SNP (to also adjust distance)
        s1 <- s2

      } else {
        # we create a new cluster with the current j SNP (s2)
        i <- j
        clustern <- names(snps)[j]
        cluster[[clustern]] <- clustern
        break
      }
    }
  }

  # for annotating the eqtl with the respective sentinel
  eqtl_subs <- cbind(eqtl_subs, sentinel=NA_character_)

  # now we select a representative SNP for each cluster using
  # 1) the number of trans associations
  # 2) the MAF
  for(i in 1:length(cluster)) {
    cl_snps <- cluster[[i]]

    # get the number of trans associations per SNP
    sub <- eqtl_subs[eqtl_subs$SNP %in% cl_snps,]

    ulength <- function(x) {return(unique(length(x)))}
    counts <- tapply(sub$GeneSymbol, sub$SNP, ulength)

    # select SNPs with maximal trans assoc
    counts <- counts[counts == max(counts)]

    # check whether we have to check MAF as well, i.e. ties in the counts
    sentinel <- NULL
    if(length(counts) > 1) {
      mafs <- MAF[names(counts)]
      # select one (random) snp with the highest MAF
      sentinel <- names(mafs[mafs == max(mafs)])[1]
    } else {
      sentinel <- names(counts)
    }

    # set the sentinel in the eQTL list
    set <- eqtl_subs$SNP %in% cl_snps
    eqtl_subs[set]$sentinel <- rep(sentinel, sum(set))
  }

  eqtl_subs
}, mc.cores = threads)

# get single data.table
full <- do.call(rbind, full)
# only sentinel SNP associations
pruned <- full[full$SNP == full$sentinel]

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=ffull, full)
saveRDS(file=fpruned, pruned)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
