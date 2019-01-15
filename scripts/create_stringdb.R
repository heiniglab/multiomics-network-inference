#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
#' Creates a rds file containing the filtered string database
#' as a graph object
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
# Load libraries and scripts
# ------------------------------------------------------------------------------
library(igraph)
library(graph)
library(data.table)

# ------------------------------------------------------------------------------
# Get snakemake parms
# ------------------------------------------------------------------------------
fout <- snakemake@output[[1]]
fstring <- snakemake@input$string
fgtex <- snakemake@input$gtex
fgene_annog <- snakemake@input$gene_annot

# ------------------------------------------------------------------------------
# Load string db as graph object
# ------------------------------------------------------------------------------
string.all <- fread(fstring,
                    data.table=F, header=T, stringsAsFactors=F)
string.inter <- string.all[string.all$experimental>=1 | string.all$database>=1,]

string.nodes <- unique(c(string.inter[,1], string.inter[,2]))
string.db <- graphNEL(nodes=string.nodes)
string.db <- addEdge(string.inter[,1], 
                     string.inter[,2], 
                     string.db)

# ------------------------------------------------------------------------------
# filtered for expressed genes in gtex whole blood
# ------------------------------------------------------------------------------
expr = read.csv(fgtex, 
                sep="\t", 
                skip=2, 
                stringsAsFactors=F)

expressed <- expr[expr[,"Whole.Blood"] > 0.1, "Name"]
ga <- load_gene_annotation(fgene_annot)
expressed.symbols <- ga[intersect(expressed, names(ga))]$SYMBOL

string.nodes = intersect(nodes(string.db), c(expressed.symbols))

string.db = subGraph(string.nodes, string.db)

# ------------------------------------------------------------------------------
# get largest connected component
# ------------------------------------------------------------------------------
ig = graph_from_graphnel(string.db)
cl = clusters(ig)
keep = nodes(string.db)[cl$membership == which.max(cl$csize)]
string.db = subGraph(keep, string.db)

# ------------------------------------------------------------------------------
# save
# ------------------------------------------------------------------------------
saveRDS(string.db, ofile)
