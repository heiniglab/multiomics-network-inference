#' Creates the prior matrix for a specific sentinel
#' 
#' @author Johann Hawe
#' 

library(qvalue)
library(data.table)
library(graph)
library(fdrtool)
library(Homo.sapiens)
library(pheatmap)

# source necessary scripts
source("scripts/lib.R")
source("scripts/priors.R")

feqtl <- snakemake@input[["gg_priors"]]
fsnpinfo <- snakemake@input[["eqtl_priors"]]
fcpgcontext <- snakemake@input[["cpg_context"]]
fexpr <- snakemake@input[["ranges"]]
fstring <- snakemake@input[["string"]]
franges <- snakemake@input[["ranges"]]
sentinel <- snakemake@params$sentinel
ofile <- snakemake@output[[1]]
fplot <- snakemake@params$plot_file

print("Loading string db.")
string_db <- readRDS(fstring)

# load ranges object
ranges <- readRDS(franges)

# get all entities as single vector
nodes <- with(ranges, c(sentinel, cpg_genes$SYMBOL, 
                        snp_genes$SYMBOL, names(cpgs)))
if(!is.null(ranges$tfs)){
  nodes <- c(nodes, ranges$tfs$SYMBOL)
}
if(!is.null(ranges$spath)){
  nodes <- c(nodes, ranges$spath$SYMBOL)
}
nodes <- unique(nodes)

# simply delegate
pr <- get_link_priors(ranges, nodes, string_db, fcpgcontext)

saveRDS(pr, file=ofile)

# create a heatmap to be able to look at the priors
pheatmap(filename = fplot, pr)
