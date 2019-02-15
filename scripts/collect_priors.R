#' Creates the prior matrix for a specific sentinel
#'
#' @author Johann Hawe
#'

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts.")
# ------------------------------------------------------------------------------
library(qvalue)
library(data.table)
library(graph)
library(fdrtool)
library(Homo.sapiens)
library(pheatmap)

source("scripts/lib.R")
source("scripts/priors.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
# input
fgene_priors <- snakemake@input[["gg_priors"]]
feqtl_priors <- snakemake@input[["eqtl_priors"]]
fcpgcontext <- snakemake@input[["cpg_context"]]
fexpr <- snakemake@input[["ranges"]]
fppi <- snakemake@input[["ppi"]]
franges <- snakemake@input[["ranges"]]
fcpg_annot <- snakemake@input[["cpg_annot"]]

# params
sentinel <- snakemake@wildcards$sentinel

# output
fout <- snakemake@output[[1]]
fplot <- snakemake@output[[2]]

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------

print("Loading PPI db.")
ppi_db <- readRDS(fppi)

# load ranges object
ranges <- readRDS(franges)

# get all entities as single vector
nodes <- c(sentinel, ranges$snp_genes$SYMBOL, 
           ranges$tfs$SYMBOL, ranges$spath$SYMBOL)

if(ranges$seed == "meqtl") {
  nodes <- c(nodes, with(ranges,
                         c(cpg_genes$SYMBOL, names(cpgs))))
} else {
  nodes <- c(nodes, ranges$trans_genes$SYMBOL)
}
nodes <- unique(nodes)

# ------------------------------------------------------------------------------
print("Load gtex priors, extract link priors.")
# ------------------------------------------------------------------------------
load.gtex.priors(sentinel, fgene_priors, feqtl_priors)

pr <- get_link_priors(ranges, nodes, ppi_db, fcpgcontext, fcpg_annot)

# create a heatmap to be able to look at the priors
pheatmap(filename = fplot, pr)

# ------------------------------------------------------------------------------
print("Saving data.")
# ------------------------------------------------------------------------------
saveRDS(pr, file=fout)

