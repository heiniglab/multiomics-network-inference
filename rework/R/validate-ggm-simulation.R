#'
#' Script to validate the GGMs calculated on simulated
#' data
#' 
#' @author Johann Hawe
#'


sink(file=snakemake@log[[1]])

# load needed libraries
library(BDgraph, lib="~/epigenereg/packages/2017/R/3.4/")
library(GenomicRanges)

source("R/lib.R")

ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]

print(paste0("Processing file: ", ifile, "."))

# load data
load(ifile)

# use the bdgraph internal method to get spec/sens and f1
perf <- t(compare(data.sim, ggm_fit, ggm_fit_no_priors))
rownames(perf) <- c("true", "ggm_fit", "ggm_fit_no_priors")

# write result to output file
write.table(file=ofile, perf, col.names=NA, sep="\t", quote=F, row.names = TRUE)

sink()