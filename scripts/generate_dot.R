#------------------------------------------------------------------------------
#' Script to generate a dot file which can be used for example
#' with Cytoscape to visualize graph models.
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
source("scripts/lib.R")

#------------------------------------------------------------------------------
print("Getting snakemake params.")
#------------------------------------------------------------------------------

# get in and output
ffit <- snakemake@input$fits
#ffit_old <- snakemake@input$old

fout <- snakemake@output[[1]]

# get wildcards
graph_type <- snakemake@wildcards$graph_type
sentinel <- snakemake@wildcards$sentinel
cohort <- snakemake@wildcards$cohort
if(is.null(cohort)) cohort <- "GTEx"

# define available graph types
gtypes <- c("bdgraph", "bdgraph_no_priors", "genenet", "irafnet", 
            "glasso", "glasso_no_priors", "genie3")

if(!graph_type %in% gtypes) {
  stop(paste0("Graph type not supported: ", graph_type))
}
print("Using graph type:")
print(graph_type)

#------------------------------------------------------------------------------
print("Loading data.")
#------------------------------------------------------------------------------
#if(graph_type %in% c("glasso", "glasso_no_priors", "genie3")) {
  fits <- readRDS(ffit)
#} else {
#  fits <- readRDS(ffit_old)
#}
g <- fits[[graph_type]]

print("Loaded graph:")
print(g)

#------------------------------------------------------------------------------
print("Creating the dot file.")
#------------------------------------------------------------------------------
attrs <- plot_ggm(g, sentinel, graph.title=paste0(c(sentinel, cohort, graph_type), 
                                                  collapse="|"), 
                  plot.on.device=F, dot.out=fout)

