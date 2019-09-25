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
source("scripts/collect_data_methods_tfa.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# input
franges <- snakemake@input[["ranges"]]
fkora_data <- snakemake@input[["kora"]]
flolipop_data <- snakemake@input[["lolipop"]]
fkora_activities <- snakemake@input$kora_activities
flolipop_activities <- snakemake@input$lolipop_activities
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
expr_syms = c(unlist(ranges$snp_genes$SYMBOL),
                unlist(ranges$cpg_genes$SYMBOL),
                unlist(ranges$tfs$SYMBOL),
                unlist(ranges$spath$SYMBOL),
                unlist(ranges$trans_genes$SYMBOL))
expr_syms <- unique(expr_syms)
meth_probes <- names(ranges$cpgs)

# ------------------------------------------------------------------------------
print("Loading cohort data.")
# ------------------------------------------------------------------------------
if("kora" %in% cohort) {
  load(fkora_data)
  acts <- readRDS(fkora_activities)
} else if ("lolipop" %in% cohort) {
  load(flolipop_data)
  acts <- readRDS(flolipop_activities)
} else {
  stop("Cohort not supported.")
}

# ------------------------------------------------------------------------------
print("Create merged data frame.")
# ------------------------------------------------------------------------------
genes <- colnames(acts)[colnames(acts) %in% expr_syms]
data <- cbind.data.frame(covars,
                         acts[,genes, drop=F],
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
print(meth_probes[!meth_probes %in% colnames(meth)])

# remove not needed remaining data frames
rm(expr, meth, covars, acts)
gc()

print(paste0("Data dimensions: ", paste(dim(data), collapse=",")))

# ------------------------------------------------------------------------------
print("Removing covariate effects from raw data.")
# ------------------------------------------------------------------------------
data <- adjust_data(sentinel, ranges, genes, data, geno, fccosmo, fceqtl)

# ------------------------------------------------------------------------------
print("Saving adjusted data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, data)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
