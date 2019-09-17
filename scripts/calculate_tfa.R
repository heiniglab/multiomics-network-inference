#' -----------------------------------------------------------------------------
#' Prepare a gene expression matrix, residual normalized for covariates, and
#' substitute TF expression values with TranscriptionFactor Activities
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Jul 31 14:18:49 2019
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(plsgenomics)
library(pheatmap)
source("scripts/lib.R")
source("scripts/collect_data_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------
# input
fcohort_data <- snakemake@input$cohort_data
ftfbs_annot <- snakemake@input$tfbs_annot

# output
fout_heatmap <- snakemake@output$heatmap
fout_tfa <- snakemake@output$tfa
fout_expr <- snakemake@output$expr

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
load(fcohort_data)

# get df containing only expr covariates and the expr probes
# needed for covariate removal method
expr_covars <- cbind.data.frame(expr,
                                covars[,c("age", "sex", "RIN",
                                          "batch1", "batch2")])
probe_resid <- rm_covariate_effects(expr_covars, "expr")

# we need one expression value per sample per gene -> summarize probes belonging
# to one gene
all_syms <- symbols.from.probeids(colnames(probe_resid))
symbol_resid <- summarize(probe_resid, all_syms)

# get the tfbs tss annotation
tss_annot <- readRDS(ftfbs_annot)
tfs <- unique(sapply(strsplit(colnames(tss_annot), "\\."), "[[", 1))

# one TF might have the same target measured more than once -> summarize
tss_annot_summarized <- sapply(tfs, function(tf) {
  rowSums(tss_annot[,grepl(paste0(tf, "\\."), colnames(tss_annot)), drop=F])
})

# ------------------------------------------------------------------------------
print("Define annotation and data subsets.")
# ------------------------------------------------------------------------------
# get TFs and their targets
targets <- intersect(rownames(tss_annot), colnames(symbol_resid))
tf_sub <- tfs[tfs %in% colnames(symbol_resid)]

# get the annotation and data subsets
annot_sub <- tss_annot_summarized[targets,,drop=F]
annot_sub <- annot_sub[, tf_sub]
data_sub <- t(symbol_resid[, targets])

# ------------------------------------------------------------------------------
print("Estimating TFAs using PLS/SIMPLS and substituting.")
# ------------------------------------------------------------------------------
TFA <- t(TFA.estimate(annot_sub, data_sub)$TFA)
colnames(TFA) <- colnames(annot_sub)
rownames(TFA) <- colnames(data_sub)

# substitute expression of TFs with TFA
symbol_resid_tfa <- symbol_resid
for(tf in tf_sub) {
  symbol_resid_tfa[,tf] <- TFA[,tf]
}

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fout_tfa, symbol_resid_tfa)
saveRDS(file=fout_expr, symbol_resid)
# ------------------------------------------------------------------------------
print("Getting correlations of TFs/Targets and plotting heatmap.")
# ------------------------------------------------------------------------------
# we look at the relevant subset to speed up computations for cor/pheatmap
gene_subset <- unique(c(targets, tf_sub))
cor_tfa <- cor(symbol_resid_tfa[, gene_subset])
cor_exp <- cor(symbol_resid[, gene_subset])

cor_comb <- cor_tfa
cor_comb[lower.tri(cor_comb, diag=F)] <- cor_exp[lower.tri(cor_exp, diag=F)]

# plot heatmaps of correlations
ph <- pheatmap(cor_comb, cluster_cols=F, cluster_rows=F)
png(fout_heatmap, res=900)
print(ph)
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")