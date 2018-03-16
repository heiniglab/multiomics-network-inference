#!/usr/bin/env Rscript

#' Creates a rds file containing the filtered string database
#' as a graph object
#'
#' @author Johann Hawe

library(igraph)
library(graph)
library(data.table)
library(Homo.sapiens)

# get needed file paths
ofile <- snakemake@output[[1]]
string_file <- snakemake@input[["string"]]
gtex_file <- snakemake@input[["gtex"]]

# load db anew
string.all <- fread(string_file,
                    data.table=F, header=T, stringsAsFactors=F)
string.inter <- string.all[string.all$experimental>=1 | string.all$database>=1,]
rm(string.all)

# create graph object
string.nodes <- unique(c(string.inter[,1], string.inter[,2]))
string.db <- graphNEL(nodes=string.nodes)
string.db <- addEdge(string.inter[,1], 
                     string.inter[,2], 
                     string.db)

# filtered for expressed genes in gtex whole blood
expr = read.csv(gtex_file, 
                sep="\t", 
                skip=2, 
                stringsAsFactors=F)

expressed = expr[expr[,"Whole.Blood"] > 0.1, "Name"]
expressed = sapply(strsplit(expressed, "." ,fixed=T) ,"[", 1)
expressed.symbols = select(Homo.sapiens, keys=expressed, keytype="ENSEMBL", columns="SYMBOL")
expressed.symbols = unique(expressed.symbols[,"SYMBOL"])

string.nodes = intersect(nodes(string.db), expressed.symbols)
string.db = subGraph(string.nodes, string.db)

# get largest connected component
ig = graph_from_graphnel(string.db)
cl = clusters(ig)
keep = nodes(string.db)[cl$membership == which.max(cl$csize)]
string.db = subGraph(keep, string.db)

# save
saveRDS(string.db, ofile)
