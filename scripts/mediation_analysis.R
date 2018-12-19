#' -----------------------------------------------------------------------------
#' In this script we evaluate mediation performance of all SNP genes for a
#' specific sentinel on the provided cohort data.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Oct 31 11:15:50 2018
#' -----------------------------------------------------------------------------
sink(snakemake@log[[1]], append = F, type = "output")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
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

# the snp genes
sgenes <- ranges$snp_genes$SYMBOL
sgenes <- sgenes[sgenes %in% colnames(data)]

# the trans associated entities
if(ranges$seed == "meqtl") {
  ta <- names(ranges$cpgs)
} else {
  ta <- ranges$trans_genes$SYMBOL
}
ta <- ta[ta %in% colnames(data)]

# ------------------------------------------------------------------------------
print("Performing mediation analysis.")
# ------------------------------------------------------------------------------
med <- mediation(data, sentinel, sgenes, ta,
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