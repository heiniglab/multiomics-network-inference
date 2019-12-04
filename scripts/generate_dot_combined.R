#------------------------------------------------------------------------------
#' Script to generate a dot file for cohort combined graphs, i.e. takes the
#' fits on the individual cohorts, merges the graphs, annotates them again
#' and then writes out the dot file
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#------------------------------------------------------------------------------
print("Loading libraries and scripts.")
#------------------------------------------------------------------------------
library(graph)
library(igraph)
source("scripts/reg_net.R")
source("scripts/lib.R")

#------------------------------------------------------------------------------
print("Getting snakemake params.")
#------------------------------------------------------------------------------

# get in and output
ffit_kora <- snakemake@input$new_kora
ffit_lolipop <- snakemake@input$new_lolipop

franges <- snakemake@input$ranges
fppi_db <- snakemake@input$ppi_db
fcpg_context <- snakemake@input$cpg_context
ftss_context <- snakemake@input$tss_context

# the output dot file and the combined graph
fout_dot <- snakemake@output$dot
fout_graph <- snakemake@output$graph

# get wildcards
graph_type <- snakemake@wildcards$graph_type
sentinel <- snakemake@wildcards$sentinel
seed <- snakemake@wildcards$seed

# use CpG context for meQTLs only
if(seed == "meqtl") {
  fcontext <- fcpg_context
} else {
  fcontext <- ftss_context
}

# define available graph types
gtypes <- c("bdgraph", "bdgraph_no_priors", "genenet", "irafnet",
            "glasso", "glasso_no_priors", "genie3")
if(!graph_type %in% gtypes) {
  stop(paste0("Graph type not supported: ", graph_type))
}
print("Using graph type:")
print(graph_type)

#-------------------------------------------------------------------------------
print("Loading data.")
#-------------------------------------------------------------------------------
ranges <- readRDS(franges)
ppi_db <- readRDS(fppi_db)

# we regenerated fits for glasso and genie3 without calculating all others
# we use 'old fits' for all other models
#if(graph_type %in% c("glasso", "glasso_no_priors", "genie3")) {
  fit_kora <- readRDS(ffit_kora)
  fit_lolipop <- readRDS(ffit_lolipop)
#} else {
#  fit_kora <- readRDS(ffit_kora_old)
#  fit_lolipop <- readRDS(ffit_lolipop_old)
#}
if(!graph_type %in% names(fit_kora) | !graph_type %in% names(fit_lolipop)) {
  stop("Graph type not available in fitting results.")
}

g_kora <- fit_kora[[graph_type]]
g_lolipop <- fit_lolipop[[graph_type]]

print("Loaded graphs:")
print("KORA:")
print(g_kora)
print("LOLIPOP:")
print(g_lolipop)

# ------------------------------------------------------------------------------
print("Combining individual cohort graphs.")
# ------------------------------------------------------------------------------
g_combined <- combine_graphs(g_kora, g_lolipop)

# filter for >0 node degrees
ds <- graph::degree(g_combined)
use <- names(ds[ds>0])
final_graph <- subGraph(use, g_combined)
final_graph_annotated <- annotate.graph(final_graph, ranges, ppi_db,
                                        fcontext = fcontext)

print("Final graph:")
print(final_graph_annotated)

#------------------------------------------------------------------------------
print("Saving output files.")
#------------------------------------------------------------------------------
attrs <- plot_ggm(final_graph_annotated, sentinel,
                  graph.title=paste0(c(sentinel, graph_type),
                                                  collapse="|"),
                  plot.on.device=F, dot.out=fout_dot)
saveRDS(file = fout_graph, final_graph_annotated)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
