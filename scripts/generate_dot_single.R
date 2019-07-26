#-------------------------------------------------------------------------------
#' Script to generate a dot file from graph fits without using snakemake
#' which can be used for example with Cytoscape to visualize graph models.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
print("Loading libraries and scripts.")
#-------------------------------------------------------------------------------
suppressPackageStartupMessages(library(graph))
source("scripts/lib.R")

#-------------------------------------------------------------------------------
print("Loading data.")
#-------------------------------------------------------------------------------
fin <- commandArgs(trailingOnly=T)[1]
gtype <- commandArgs(trailingOnly=T)[2]
print(paste0("Graph-type is: ", gtype))

#-------------------------------------------------------------------------------
print("Loading data.")
#-------------------------------------------------------------------------------
fits <- readRDS(fin)
g <- fits[[gtype]]

print("Graph:")
print(g)

#-------------------------------------------------------------------------------
print("Creating dot file.")
#-------------------------------------------------------------------------------
fout <- paste0(gtype, ".dot")
attrs <- plot_ggm(g, plot.on.device=F, dot.out=fout)

