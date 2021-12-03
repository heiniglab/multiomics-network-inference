#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
#' Apply graphical model inference on data collected for a specific sentinel.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Prepare libraries and source scripts.")
# ------------------------------------------------------------------------------
library(pheatmap)
suppressPackageStartupMessages(library(GenomicRanges))
library(igraph)
library(graph)
library(reshape2)
source("scripts/lib.R")
source("scripts/reg_net.R")
source("scripts/reg_net_utils.R")

# ------------------------------------------------------------------------------
print("Get snakemake parameters.")
# ------------------------------------------------------------------------------

#input
franges <- snakemake@input[["ranges"]]
fdata <- snakemake@input[["data"]]
fppi_db <- snakemake@input[["ppi_db"]]
fpriors <- snakemake@input[["priors"]]
fcpg_context <- snakemake@input[["cpg_context"]]
ftss_context <- snakemake@input[["tss_context"]]

# output
fout <- snakemake@output$fit
fsummary_plot <- snakemake@output$summary_file

# params
threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Load and prepare data.")
# ------------------------------------------------------------------------------
data <- readRDS(fdata)

# remove (rare) all-NA cases. This can happen due to scaling of all-zero entities,
# which can arise due to a very large number of cis-meQTLs which effects get
# removed from the CpGs during data preprocessing.
# NOTE: we could possibly handle this differently? Seems that these particular
# cpgs are highly influenced by genetic effects?
use <- apply(data,2,function(x) (sum(is.na(x)) / length(x)) < 1)
data <- data[,use]

print("Dimensions of data:")
print(dim(data))

priors <- readRDS(fpriors)

# filter for available data in priros
priors <- priors[colnames(data), colnames(data)]
ranges <- readRDS(franges)

# load PPI DB
ppi_db <- readRDS(fppi_db)

# ------------------------------------------------------------------------------
# for catching the genenet summary plots, we open the pdf connection here
# ------------------------------------------------------------------------------
pdf(fsummary_plot)

# ------------------------------------------------------------------------------
print("Infer regulatory networks.")
# ------------------------------------------------------------------------------
if(ranges$seed == "meqtl") {
  fcontext <- fcpg_context
} else {
  fcontext <- ftss_context
}

result <- infer_all_graphs(data, priors, ranges, fcontext, ppi_db,
                           threads)

dev.off()

# ------------------------------------------------------------------------------
print("All done. Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, result)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
