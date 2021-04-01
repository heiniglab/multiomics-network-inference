#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
#' Apply graphical model inference on data collected for a specific sentinel.
#' Adds a certain amount of noise to the prior information prior to inference.
#'
#' @author Johann Hawe <johann.hawe@tum.de>
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

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

# ------------------------------------------------------------------------------
print("Get snakemake parameters.")
# ------------------------------------------------------------------------------

#input
franges <- snakemake@input$ranges
fdata_kora <- snakemake@input$data_kora
fdata_lolipop <- snakemake@input$data_lolipop
fppi_db <- snakemake@input$ppi_db
fpriors <- snakemake@input$priors
fcpg_context <- snakemake@input$cpg_context
ftss_context <- snakemake@input$tss_context

# output
fout <- snakemake@output$fit
fsummary_plot <- snakemake@output$summary_file

# params
threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Load and prepare data.")
# ------------------------------------------------------------------------------
remove_all_na <- function(data) {
  # remove (rare) all-NA cases. This can happen due to scaling of all-zero entities,
  # which can arise due to a very large number of cis-meQTLs which effects get
  # removed from the CpGs during data preprocessing.
  # NOTE: we could possibly handle this differently? Seems that these particular
  # cpgs are highly influenced by genetic effects?
  use <- apply(data, 2, function(x)
    (sum(is.na(x)) / length(x)) < 1)
  data <- data[, use]
  data
}

data_kora <- remove_all_na(readRDS(fdata_kora))
data_lolipop <- remove_all_na(readRDS(fdata_lolipop))

print("Dimensions of KORA data:")
print(dim(data_kora))

print("Dimensions of LOLIPOP data:")
print(dim(data_lolipop))

# we only look at replication, so we only consider the nodes present in both
# cohorts
common_nodes <- intersect(colnames(data_kora), colnames(data_lolipop))


# filter for available data in priros
priors <- readRDS(fpriors)
priors <- priors[common_nodes, common_nodes]
ranges <- readRDS(franges)

# load PPI DB
ppi_db <- readRDS(fppi_db)

if (ranges$seed == "meqtl") {
  fcontext <- fcpg_context
} else {
  fcontext <- ftss_context
}

# ------------------------------------------------------------------------------
# for catching the genenet summary plots, we open the pdf connection here
# ------------------------------------------------------------------------------
pdf(fsummary_plot)

# ------------------------------------------------------------------------------
# We add different levels of noise and infer graphs for all scenarios
print("Infer regulatory networks.")
# ------------------------------------------------------------------------------

# helper to get a noisy prior matrix
noisify_priors <- function(priors, noise_level) {
  
  if(noise_level == 0) return(priors)
  
  PSEUDO_PRIOR <- 1e-7
  number_edges_with_prior <- sum(priors[upper.tri(priors)] > PSEUDO_PRIOR)
  number_edges_without_prior <- sum(priors[upper.tri(priors)] == PSEUDO_PRIOR)
  
  number_entries_to_switch <-
    round(noise_level * number_edges_with_prior)
  
  # rare cases where we have more prior annotated edges than no-prior edges
  if(number_entries_to_switch > number_edges_without_prior) {
    number_entries_to_switch <- number_edges_without_prior
  }
  
  prior_idx <-
    sample(which(upper.tri(priors) & priors > PSEUDO_PRIOR),
           number_entries_to_switch)
  non_prior_idx <-
    sample(which(upper.tri(priors) & priors == PSEUDO_PRIOR),
           min(number_entries_to_switch, ))
  
  # swap idxs
  temp <- priors[prior_idx]
  priors[prior_idx] <- priors[non_prior_idx]
  priors[non_prior_idx] <- temp
  
  return(priors)
}

# include 0 noise ('normal' model)
noise_levels <- seq(0, 0.8, by = 0.2)
result <- lapply(noise_levels, function(noise_level) {
  priors <- noisify_priors(priors, noise_level)
  result_kora <-
    infer_all_graphs(data_kora, priors, ranges, fcontext, ppi_db,
                            threads)
  result_lolipop <-
    infer_all_graphs(data_lolipop, priors, ranges, fcontext, ppi_db,
                            threads)
  list(kora = result_kora, lolipop = result_lolipop)
})
names(result) <- paste0("noise_level_", noise_levels)
dev.off()

# ------------------------------------------------------------------------------
print("All done. Saving results.")
# ------------------------------------------------------------------------------
saveRDS(file = fout, result)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
