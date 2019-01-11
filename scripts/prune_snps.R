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
  snps <- with(eqtl_subs, GRanges(paste0("chr", chr,
                                         IRanges(SNPPos, width=1))))
  names(snps) <- eqtl_subs$SNP
  snps <- unique(snps)

  # get KORA genotypes
  geno <- scan_snps(snps, fdosage, individuals)

  # TODO cope with missing genotypes for some SNPs

  # calculate pairwise r^2
  pw_r2 <- cor(t(geno))^2

  # check each SNP pair (within distance threshold) for large r^2
  for(i in 1:length(snps)) {
    for (j in (i + 1):length(snps)) {
      s1 <- snps[i]
      s2 <- snps[j]
      r2 <-pw_r2[names(s1), names(s2)]

      #
      if (r2 > 0.2 & distance(s1,s2) < ld_max_dist) {

      }
    }
  }
})


# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
