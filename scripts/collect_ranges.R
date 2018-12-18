#'
#' Script to collect all associated ids (expr, meth and geno)
#' for a specific sentinel snp
#'
#' @author Johann Hawe
#'
#' @date 2017/03/04
#'

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
# load packages, source scripts
# ------------------------------------------------------------------------------
library(GenomicRanges)
library(GenomicFeatures)
library(FDb.InfiniumMethylation.hg19)
library(data.table)
library(illuminaHumanv3.db)
library(rtracklayer)
library(graph)
library(RBGL) # for shortest paths
library(Matrix)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# Define needed methods
# ------------------------------------------------------------------------------

#' Gets the shortest paths between two sets of genes
#'
#' Uses the validated string network and identified the genes on the shortest paths
#' between two given genesets.
#'
#' @param cis List of nodes in cis (e.g. CpGs)
#' @param trans List of nodes in trans (e.g. snp genes and TFs)
#' @param snp_genes List of snp/trans genes near the sentinel
#' @param best_trans The identified 'best' trans genes in the priorization
#' analysis
#' @param string_db Instance of the string.db to be used (graphNEL)
#' @return A vector of gene symbols being on the shortest path between the two lists
#' of genes as found in the validated string network
#'
#' @author Johann Hawe, Matthias Heinig
#'
get_string_shortest_paths <- function(cis, trans, snp_genes,
                                      best_trans, string_db) {

  print("Preprocessing.")
  g <- string_db
  g.nodes <- nodes(string_db)

  # ensure to have only nodes in our giant cluster
  cis <- cis[which(cis %in% g.nodes)]
  trans <- trans[which(trans %in% g.nodes)]
  snp_genes <- snp_genes[which(snp_genes %in% g.nodes)]
  if(length(cis) == 0 | length(trans) == 0 | length(snp_genes) == 0) {
    return(NULL)
  }

  # calculate weights for nodes
  prop = propagation(graph2sparseMatrix(g), n.eigs=500,
                     from=cis, to=trans, sum="both")

  if(is.null(best_trans)) {
    warning("No best trans genes detected, using propagation results.")
    # get the best trans gene
    best_trans = snp_genes[which.max(prop[snp_genes,"from"])]
  }

  print(paste0("Best trans: ", paste(best_trans, collapse = ",")))

  ## find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  node.weights <- rowSums(prop)
  node.weights = max(node.weights) - node.weights + 1

  print("Getting minimal node weight path.")
  # extract shortest paths
  sp = min.node.weight.path(g, node.weights, from=cis, to=best_trans)
  nodes <- setdiff(unlist(lapply(sp, "[", "path_detail")), NA)
  nodes <- setdiff(nodes, c(trans, cis))

  print("Shortest path genes:")
  print(nodes)

  return(nodes)
}

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------
fcosmo <- snakemake@input[["tcosmo"]]
fmeqtl <- snakemake@input[["meqtl"]]
fstring <- snakemake@input[["string"]]
fprio_tab <- snakemake@input$priorization

# TODO: create this file from scratch!
fcpgcontext <- snakemake@input[["cpgcontext"]]

ofile <- snakemake@output[[1]]
sentinel <- snakemake@wildcards$sentinel

# ------------------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------------------

print("Loading data.")

gene_annot <- get.gene.annotation()
gene_annot$ids <- probes.from.symbols(gene_annot$SYMBOL,
                                           as.list=T)

string_db <- readRDS(fstring)

# load trans-meQTL table
trans_meQTL = read.csv(fmeqtl,
                       sep="\t", stringsAsFactors=F)

# load trans-cosmo information
cosmo <- readRDS(fcosmo)


# load priorization table
prio <- read.table(fprio_tab, sep="\t", header=T,stringsAsFactors = F)

# get trans-genes which should be used for shortest path extraction
prio <- prio[prio$sentinel == sentinel,,drop=F]
if(nrow(prio) > 0) {
  best_trans <- unique(prio$trans.gene)
} else {
  best_trans <- NULL
}

# ------------------------------------------------------------------------------
# Collect and save ranges
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Collecting SNP and CpG ranges.")
# ------------------------------------------------------------------------------

pairs = which(trans_meQTL[,"sentinel.snp"] == sentinel)

# get sentinel idx
sentinel_idx <- which(cosmo$snp==sentinel)[1]

# the large interval range for the sentinel
chr <- paste0("chr", trans_meQTL[pairs,"chr.snp"][1])
start <- trans_meQTL[pairs,"interval.start.snp"][1]
end <- trans_meQTL[pairs,"interval.end.snp"][1]

sentinel_extrange <- GRanges(chr, IRanges(start,end))
sentinel_range <- with(cosmo[sentinel_idx,],
                       GRanges(paste0("chr", snp.chr),
                               IRanges(snp.pos, width=1)))
names(sentinel_range) <- names(sentinel_extrange) <- sentinel

# get cosmo subset
idxs <- get.trans.cpgs(sentinel, trans_meQTL, cosmo)

# get related genes, i.e. genes near meQTL loci (snp+cpg)

# extend cpgs
cosmosub <- cosmo[idxs,]

croi <- with(cosmosub, GRanges(paste0("chr", cpg.chr),
                               IRanges(cpg.pos,width=2)))
names(croi) <- as.character(cosmosub[,"cpg"])
croi <- unique(croi)

# extended sentinel region
sroi <- sentinel_extrange
names(sroi) <- sentinel

# ------------------------------------------------------------------------------
print("Retrieving SNP and CpG genes.")
# ------------------------------------------------------------------------------

# get the relevant snp genes by overlapping with our sentinel region
genes_sroi <- subsetByOverlaps(gene_annot, sroi, ignore.strand=T)

# get genes near our cpg regions
genes_by_cpg <- get.nearby.ranges(croi, promoters(gene_annot))
# get original ranges (not promoters)
genes_by_cpg <- lapply(names(genes_by_cpg), function(cg) {
  gs <- genes_by_cpg[[cg]]
  gs <- gene_annot[gene_annot$SYMBOL %in% gs$SYMBOL]
  gs
})
names(genes_by_cpg) <- names(croi)

# get single list of all cpg genes
genes_croi <- unique(unlist(GRangesList(unlist(genes_by_cpg))))

# ------------------------------------------------------------------------------
print("Collecting TFs and shortest path genes.")
# ------------------------------------------------------------------------------

tfs <- NULL
sp <- NULL

# get cpg ids and SNP gene symbols
cpgs <- names(croi)
snp_genes <- unique(genes_sroi$SYMBOL)

# modify string_db to contain our CpGs

# load the cpg-tf context
tfbs_ann <- get.chipseq.context(names(croi), fcpgcontext)
cpgs_with_tfbs <- cpgs[cpgs %in% rownames(tfbs_ann[rowSums(tfbs_ann)>0,])]
snp_genes_in_string <- snp_genes[snp_genes %in% nodes(string_db)]

# get locus graph
locus_graph <- add.to.graphs(list(string_db), sentinel, snp_genes,
                             cpgs_with_tfbs, tfbs_ann)[[1]]

# get tfs connected to cpgs
tfs = unique(unlist(adj(locus_graph, cpgs_with_tfbs)))
print(paste0("Annotated TFs: ", paste(tfs, collapse=", ")))

if(length(tfs)<1){
  warning("No TFs, skipping shortest paths calculation.")
} else {
  # the nodes we want to keep
  nodeset <- c(nodes(string_db), setdiff(tfs, "KAP1"),
               snp_genes_in_string, cpgs_with_tfbs)
  locus_graph <- subGraph(intersect(nodes(locus_graph), nodeset), locus_graph)

  syms_sp <- get_string_shortest_paths(cis = cpgs_with_tfbs,
                                       trans=unique(c(snp_genes_in_string,
                                                      tfs)),
                                       snp_genes=snp_genes_in_string,
                                       best_trans,
                                       locus_graph)

  if(length(syms_sp) < 1){
    warning("No shortest path genes.")
  } else {
    sp <- gene_annot[gene_annot$SYMBOL %in% syms_sp]
    sp$ids <- probes.from.symbols(sp$SYMBOL, as.list=T)
  }
  tfs <- gene_annot[gene_annot$SYMBOL %in% tfs]
  tfs$ids <- probes.from.symbols(tfs$SYMBOL, as.list=T)
}

# ------------------------------------------------------------------------------
print("Finalizing and saving results.")
# ------------------------------------------------------------------------------
result <-  list(cpgs=croi,sentinel=sentinel_range,
                sentinel_ext_range=sentinel_extrange,
                snp_genes=genes_sroi, cpg_genes=genes_croi,
                cpg_genes_by_cpg=genes_by_cpg)

if(!is.null(sp)){
  result$spath <- sp
}
if(!is.null(tfs)){
  result$tfs <- tfs
}
# set seed name
result$seed <- "meqtl"

saveRDS(file=ofile, result)
