#' -----------------------------------------------------------------------------
#' Create a graph object reflecting BIOGRID PPIs which can be used for
#' further analysis.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Nov 12 16:22:00 2018
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(igraph))
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# ppis
fbiogrid <- snakemake@input$biogrid
fgene_annot <- snakemake@input$gene_annot

# gtex gene expression data (we filter genes for the ones expressed in blood)
fgtex <- snakemake@input$gtex

fout <- snakemake@output[[1]]

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
biogrid <- fread(fbiogrid, stringsAsFactors=F)
gtex <- fread(fgtex, stringsAsFactors=F)

# ------------------------------------------------------------------------------
print("Creating interaction network.")
# ------------------------------------------------------------------------------
nodes <- unique(c(biogrid$`Official Symbol Interactor A`,
                biogrid$`Official Symbol Interactor B`))
g <- graphNEL(nodes=nodes)
g <- addEdge(biogrid$`Official Symbol Interactor A`,
                     biogrid$`Official Symbol Interactor B`,
                     g)

# ------------------------------------------------------------------------------
print("Filtering PPI for expressed genes.")
# ------------------------------------------------------------------------------
expressed <- unlist(gtex[`Whole Blood` > 0.1,"Name"])
ga <- load_gene_annotation(fgene_annot)
expressed.symbols <- ga[intersect(expressed, names(ga))]$SYMBOL

nodes_to_keep = intersect(nodes(g), c(expressed.symbols))
g = subGraph(nodes_to_keep, g)

# ------------------------------------------------------------------------------
print("Getting largest connected component.")
# ------------------------------------------------------------------------------
ig = graph_from_graphnel(g)
cl = clusters(ig)
keep = nodes(g)[cl$membership == which.max(cl$csize)]
g = subGraph(keep, g)

print("Largest CC in BIOGRID filtered for blood expressed genes:")
print(g)

# ------------------------------------------------------------------------------
print("All done. Saving output.")
# ------------------------------------------------------------------------------
saveRDS(g, fout)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
