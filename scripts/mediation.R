#' -----------------------------------------------------------------------------
#' In this script we evaluate mediation performance of all SNP genes for a
#' specific sentinel on the provided cohort data.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Oct 31 11:15:50 2018
#' -----------------------------------------------------------------------------
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")
source("scripts/mediation_methods.R")

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
# remove (rare) all-NA cases. This can happen due to scaling of all-zero entities,
# which can arise due to a very large number of cis-meQTLs which effects get
# removed from the CpGs during data preprocessing.
# NOTE: we could possibly handle this differently? Seems that these particular
# cpgs are highly influenced by genetic effects?
use <- apply(data,2,function(x) (sum(is.na(x)) / length(x)) < 1)
data <- data[,use]

ranges <- readRDS(franges)

# the snp genes
sgenes <- ranges$snp_genes$SYMBOL
sgenes <- sgenes[sgenes %in% colnames(data)]
if(length(sgenes) == 0) {
  warning("No SNP genes in data matrix.")
  saveRDS(file=fout, NULL)
  q()
}
# the trans associated entities
if(ranges$seed == "meqtl") {
  ta <- names(ranges$cpgs)
} else {
  ta <- ranges$trans_genes$SYMBOL
}
ta <- ta[ta %in% colnames(data)]
print("Number of trans entities:")
print(length(ta))
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
sink(type="message")
