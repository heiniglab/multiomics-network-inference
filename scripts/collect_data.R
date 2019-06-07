# ------------------------------------------------------------------------------
#' Collects methylation, expression and genotype data for a either KORA or
#' LOLIPOP cohort based on a ranges object
#'
#' @autho Johann Hawe
#'
#' @date 02/06/2017
#'
#' @export
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


# ------------------------------------------------------------------------------
print("Load libraries, source scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")
source("scripts/collect_data_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# input
franges <- snakemake@input[["ranges"]]
fkora_data <- snakemake@input[["kora"]]
flolipop_data <- snakemake@input[["lolipop"]]
fccosmo <- snakemake@input[["ccosmo"]]
fceqtl <- snakemake@input[["ceqtl"]]

# params
cohort <- snakemake@wildcards$cohort
sentinel <- snakemake@wildcards$sentinel
seed <- snakemake@wildcards$seed

# output
fout <- snakemake@output[[1]]
fout_raw <- snakemake@output[[2]]

print(paste0("Sentinel is ", sentinel, "."))

# ------------------------------------------------------------------------------
print("Prepare probe and snp ids.")
# ------------------------------------------------------------------------------
ranges <- readRDS(franges)
expr_probes = c(unlist(ranges$snp_genes$ids),
                unlist(ranges$cpg_genes$ids),
                unlist(ranges$tfs$ids),
                unlist(ranges$spath$ids),
                unlist(ranges$trans_genes$ids))
expr_probes <- unique(expr_probes)
meth_probes <- names(ranges$cpgs)

# ------------------------------------------------------------------------------
print("Loading cohort data.")
# ------------------------------------------------------------------------------
if("kora" %in% cohort) {
  load(fkora_data)
} else if ("lolipop" %in% cohort) {
  load(flolipop_data)
} else {
  stop("Cohort not supported.")
}

# ------------------------------------------------------------------------------
print("Create merged data frame.")
# ------------------------------------------------------------------------------
data <- cbind.data.frame(covars,
                         expr[,colnames(expr) %in% expr_probes, drop=F],
                         meth[,colnames(meth) %in% meth_probes, drop=F],
                         geno[,colnames(geno) %in% sentinel, drop=F],
                         geno_ids=rownames(geno),
                         stringsAsFactors=F)

# ------------------------------------------------------------------------------
print("Save raw data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout_raw, data)

# check which probes we didnt have available
print("Unavailable probes:")
print(expr_probes[!expr_probes %in% colnames(expr)])
print(meth_probes[!meth_probes %in% colnames(meth)])

# remove not needed remaining data frames
rm(expr, meth, covars)
gc()

print(paste0("Data dimensions: ", paste(dim(data), collapse=",")))

# ------------------------------------------------------------------------------
print("Removing covariate effects from raw data.")
# ------------------------------------------------------------------------------
data <- adjust_data(sentinel, ranges, data, geno, fccosmo, fceqtl)

# ------------------------------------------------------------------------------
print("Saving adjusted data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, data)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
