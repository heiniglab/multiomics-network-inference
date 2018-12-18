#' -----------------------------------------------------------------------------
#' In this script we evaluate mediation performance of all SNP genes for a
#' specific sentinel on the provided cohort data.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Oct 31 11:15:50 2018
#' -----------------------------------------------------------------------------
sink(snakemake@log[[1]], append = F, type = "output")

#######
## TODO adjust this script to handle eQTL seeded ranges (mediation w.r.t. trans
# eQTLs rather than meQTLs)
#######

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(GenomicRanges)
source("scripts/lib.R")
source("scripts/validation.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------
# inputs
fdata <- snakemake@input$data
franges <- snakemake@input$ranges

# params
sentinel <- snakemake@wildcards$sentinel

# outputs
fbetas_per_gene_plot <- snakemake@output$betas_per_gene
fbeta_table <- snakemake@output$beta_table
fout <- snakemake@output$mediation

# ------------------------------------------------------------------------------
print("Loading and preparing data.")
# ------------------------------------------------------------------------------
data <- readRDS(fdata)
ranges <- readRDS(franges)

sgenes <- ranges$snp_genes$SYMBOL
sgenes <- sgenes[sgenes %in% colnames(data)]
cpgs <- names(ranges$cpgs)
cpgs <- cpgs[cpgs %in% colnames(data)]

# ------------------------------------------------------------------------------
print("Performing mediation analysis.")
# ------------------------------------------------------------------------------
med <- mediation(data, sentinel, sgenes, cpgs,
                  fbeta_table, fbetas_per_gene_plot)

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, med)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

sink()