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
library(ggplot2)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())

source("scripts/lib.R")
source("scripts/collect_data_methods.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------
# input
fcohort_data <- snakemake@input$cohort_data
ftfbs_annot <- snakemake@input$tfbs_annot

# output
fout_plot <- snakemake@output$plot
fout_tfa <- snakemake@output$tfa
fout_expr <- snakemake@output$expr

cohort <- snakemake@wildcards$cohort

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
symbol_resid_tfa[,tf_sub] <- TFA[,tf_sub]

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file=fout_tfa, symbol_resid_tfa)
saveRDS(file=fout_expr, symbol_resid)

# ------------------------------------------------------------------------------
print("Getting correlations of TFs/Targets and plotting.")
# ------------------------------------------------------------------------------

# get plotting data frame
data <- sapply(tf_sub, function(g) { 
  cor(symbol_resid_tfa[,g], symbol_resid[,g])
})

toplot <- data.frame(correlation=unlist(data))

pdf(fout_plot)
ggplot(toplot, aes(x=correlation)) +
  geom_histogram() +
  labs(title=paste0("Correlation between TFA and Expression in ", cohort))
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")