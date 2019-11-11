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
tfa_or_expr <- snakemake@params$tfa_or_expr

print(paste0("Using '", tfa_or_expr, "' for gene measurements"))
      
# we source a different methods script in case we get TFA
if("tfa" %in% tfa_or_expr) {
  source("scripts/collect_data_methods_tfa.R")
} else {
  source("scripts/collect_data_methods.R")
}

# output
fout <- snakemake@output[[1]]
fout_raw <- snakemake@output[[2]]

print(paste0("Sentinel is ", sentinel, "."))

# ------------------------------------------------------------------------------
print("Prepare probe and snp ids.")
# ------------------------------------------------------------------------------
ranges <- readRDS(franges)

# get expression probes
expr_probes = c(unlist(ranges$snp_genes$ids),
                unlist(ranges$cpg_genes$ids),
                unlist(ranges$tfs$ids),
                unlist(ranges$spath$ids),
                unlist(ranges$trans_genes$ids))
expr_probes <- unique(expr_probes)

# get expr gene symbols (for TFA based data)
expr_syms = c(unlist(ranges$snp_genes$SYMBOL),
                unlist(ranges$cpg_genes$SYMBOL),
                unlist(ranges$tfs$SYMBOL),
                unlist(ranges$spath$SYMBOL),
                unlist(ranges$trans_genes$SYMBOL))
expr_syms <- unique(expr_syms)

# get methylation probes
meth_probes <- names(ranges$cpgs)

# ------------------------------------------------------------------------------
print("Loading cohort data.")
# ------------------------------------------------------------------------------
if("kora" %in% cohort) {
  load(fkora_data)
  if("tfa" %in% tfa_or_expr) {
    acts <- readRDS(fkora_activities)
  }
} else if ("lolipop" %in% cohort) {
  load(flolipop_data)
  if("tfa" %in% tfa_or_expr) {
    acts <- readRDS(flolipop_activities)
  }
} else {
  stop("Cohort not supported.")
}

# ------------------------------------------------------------------------------
print("Create merged data frame.")
# ------------------------------------------------------------------------------
if("tfa" %in% tfa_or_expr) {
  print("Collecting TFA based data.")
  genes <- colnames(acts)[colnames(acts) %in% expr_syms]
  data <- cbind.data.frame(covars,
                           acts[,genes, drop=F],
                           meth[,colnames(meth) %in% meth_probes, drop=F],
                           geno[,colnames(geno) %in% sentinel, drop=F],
                           geno_ids=rownames(geno),
                           stringsAsFactors=F)
} else {
  print("Collecting EXPR based data.")
  data <- cbind.data.frame(covars,
                           expr[,colnames(expr) %in% expr_probes, drop=F],
                           meth[,colnames(meth) %in% meth_probes, drop=F],
                           geno[,colnames(geno) %in% sentinel, drop=F],
                           geno_ids=rownames(geno),
                           stringsAsFactors=F)
}

# ------------------------------------------------------------------------------
print("Saving raw data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout_raw, data)

# remove not needed remaining data frames
rm(expr, meth, covars)
gc()

print(paste0("Data dimensions: ", paste(dim(data), collapse=",")))

# ------------------------------------------------------------------------------
print("Removing covariate effects from raw data.")
# ------------------------------------------------------------------------------
if("tfa" %in% tfa_or_expr) {
  data <- adjust_data(sentinel, ranges, genes, data, geno, fccosmo, fceqtl)  
} else {
  data <- adjust_data(sentinel, ranges, data, geno, fccosmo, fceqtl)
}

# ------------------------------------------------------------------------------
print("Saving adjusted data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, data)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
