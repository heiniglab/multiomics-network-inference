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
ffit <- snakemake@input[[1]]
fout <- snakemake@output[[1]]

# get wildcards
graph_type <- snakemake@wildcards$graph_type
sentinel <- snakemake@wildcards$sentinel

# define available graph types
gtypes <- c("graph", "graph_no_priors", "genenet")
if(!graph_type %in% gtypes) {
  stop(paste0("Graph type not supported: ", graph_type))
}

#------------------------------------------------------------------------------
print("Loading data.")
#------------------------------------------------------------------------------
fits <- readRDS(ffit)
g <- fits[[graph_type]]

print("Loaded graph:")
print(g)

#------------------------------------------------------------------------------
print("Creating the dot file.")
#------------------------------------------------------------------------------
attrs <- plot.ggm(g, sentinel, plot.on.device=F, dot.out=fout)

