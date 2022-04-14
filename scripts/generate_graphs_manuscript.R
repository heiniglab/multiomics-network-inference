#' -----------------------------------------------------------------------------
#' Generate the graphs 'replicated' from the meQTL study
#' to be layouted and put into the GGM manuscript
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Nov 26 14:38:26 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(graph)
library(tidyverse)
library(Rgraphviz)
source("scripts/lib.R")
source("scripts/generate_graphs_manuscript_methods.R")

# ------------------------------------------------------------------------------
# Plot replicated graph (meQTL study vs current study)
loci <- c("rs730775")

# load our graphs
ggm <- sapply(loci, function(l) {
  parse_dot(paste0("results/current/biogrid_stringent/graph_plots_tfa/",
                   l,
                   "_meqtl/glasso_combined.dot"))
})

# load the ones from the meQTL paper
rw <- sapply(loci, function(l) {
  parse_dot(paste0("../meQTLs/results/current/rw_ggm_integration/",
                   l,
                   "/graph-final.dot"))
})

# load ranges information so that we can annotate the gene types
ranges <- sapply(loci, function(l) {
  readRDS(paste0("results/current/biogrid_stringent/ranges/",
                 l,
                 "_meqtl.rds"))
})

merged <- lapply(loci, function(l) {
 lggm <- ggm[[l]]
 lrw <- rw[[l]]
 
 g <- merge_graph(g1, g2)
 g
})

plot_graph(merged[[1]], "rs730775_tfa.pdf", "rs730775_tfa.dot", ranges[[1]])


# ------------------------------------------------------------------------------
# Comparison between graphs obtained from HURI PPI and our Biogrid PPI

locus <- "rs9274623_eqtlgen"
huri <-
  parse_dot(
    paste0(
      "results/current/huri/graph_plots_tfa/_rerun/",
      locus,
      "/glasso_combined.dot"
    )
  )
biogrid <-
  parse_dot(
    paste0(
      "results/current/biogrid_stringent/graph_plots_tfa/_rerun/",
      locus,
      "/glasso_combined.dot"
    )
  )
biogrid_ranges <-
  readRDS(paste0("results/current/biogrid_stringent/ranges/",
                 locus,
                 ".rds"))

merged <- merge_graph(huri, biogrid)

plot_graph(merged, 
           "results/current/revisions/graph_comparison_huri_biogrid.pdf", 
           "results/current/revisions/graph_comparison_huri_biogrid.dot", 
           biogrid_ranges)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
