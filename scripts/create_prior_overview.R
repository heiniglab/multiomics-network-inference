#' -----------------------------------------------------------------------------
#' Create a few overview plots for the generated priors.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Sep 20 11:18:36 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())

# TODO: use snakemake paths
gg <- readRDS("results/current/biogrid/gtex.gg.cors.rds")
eqtl <- readRDS("results/current/biogrid/gtex.eqtl.priors.rds")
tfbs <- readRDS("results/current/tfbs_tss_annot.rds")
tfs <- unique(gsub("\\..*", "", colnames(tfbs)))
gg_tf_sub <- subset(gg, g1 %in% tfs & g2 %in% tfs)

# show some general information
print(quantile(gg$prior))
print(quantile(gg_tf_sub$prior))

# plotting
toplot <- melt(gg[,c("pval", "prior")])
ggplot(toplot, aes(x=value, color=variable)) + geom_freqpoly()
       
# eqtls
eqtl$prior <- (1 - eqtl$lfdr)
toplot <- melt(eqtl[,c("pval_nominal", "prior")])
ggplot(toplot, aes(x=value, color=variable)) + geom_freqpoly()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
