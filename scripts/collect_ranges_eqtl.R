# ------------------------------------------------------------------------------
#'
#' Script to collect all associated ids (expr and geno)
#' for a specific sentinel snp of the eQTL Gen results
#'
#' @author Johann Hawe
#'
#' @date 2017/03/04
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load packages, source scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
library(data.table)
library(graph)
source("scripts/lib.R")
source("scripts/biomaRt.R")

# ------------------------------------------------------------------------------
# Define needed methods
# ------------------------------------------------------------------------------

#' Gets the shortest paths between two sets of genes
#'
#' Uses the validated string network and identified the
#' genes on the shortest paths between two given genesets.
#'
#' @param cis List of nodes in cis
#' @param trans List of nodes in trans
#' @param snp_genes List of snp genes
#' @param ppi_db Instance of the PPI database to be used (graphNEL)
#'
#' @return A vector of gene symbols being on the shortest path between the
#' two lists of genes as found in the validated string network
#'
#' @author Johann Hawe, Matthias Heinig
#'
get_shortest_paths <- function(cis, trans, snp_genes, ppi_db) {

  ppi_genes <- nodes(ppi_db)
  # ensure to have only nodes in our giant cluster
  cis <- cis[which(cis %in% ppi_genes)]
  trans <- trans[which(trans %in% ppi_genes)]
  snp_genes <- snp_genes[which(snp_genes %in% ppi_genes)]
  if(length(cis) == 0 | length(trans) == 0 | length(snp_genes) == 0) {
    return(NULL)
  }

  # calculate weights for nodes
  prop = propagation(graph2sparseMatrix(ppi_db), n.eigs=500,
                     from=cis, to=trans, sum="both")

  # get the best snp gene
  best_snp_gene = snp_genes[which.max(prop[snp_genes,"from"])]

  print(paste0("Best cis: ", paste(best_snp_gene, collapse = ",")))

  ## find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  node.weights <- rowSums(prop)
  node.weights = max(node.weights) - node.weights + 1

  print("Getting minimal node weight path.")
  # extract shortest paths
  sp = min.node.weight.path(ppi_db, node.weights, from=trans, to=best_snp_gene)
  nodes <- setdiff(unlist(lapply(sp, "[", "path_detail")), NA)
  nodes <- setdiff(nodes, c(trans, cis))

  print("Shortest path genes:")
  print(nodes)

  return(nodes)
}

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
# inputs
feqtl <- snakemake@input$eqtl
fppi <- snakemake@input$ppi
ftfbs_annot <- snakemake@input$tfbs_annot

# output file
fout <- snakemake@output[[1]]

# params
sentinel <- snakemake@wildcards$sentinel

# ------------------------------------------------------------------------------
print("Load and preprocess data.")
# ------------------------------------------------------------------------------

# PPIs
ppi_db <- readRDS(fppi)
ppi_genes <- nodes(ppi_db)

# gene annotation
gene_annot <- get.gene.annotation()
gene_annot$ids <- probes.from.symbols(gene_annot$SYMBOL,
                                      as.list=T)

# load trans-eQTL
eqtl = fread(paste0("zcat ", feqtl))
eqtl <- eqtl[eqtl$SNP == sentinel,]
if(nrow(eqtl)<5) stop("Will only process hotspots with at least 5 associations.")

# load TFBS annotation
tfbs <- readRDS(ftfbs_annot)

# ------------------------------------------------------------------------------
print("Collecting SNP and gene ranges.")
# ------------------------------------------------------------------------------

# get sentinel region + extended region in which to look for genes
spos <- get_snpPos_biomart(sentinel)
if(nrow(spos)<1) stop("Couldn't get SNP pos from biomaRt.")

sentinel_range <- with(spos,
                       GRanges(paste0("chr", chr), IRanges(start, width=1),
                               id=snp))
sentinel_extrange <- resize(sentinel_range, 1e6, fix="center")
names(sentinel_range) <- names(sentinel_extrange) <- sentinel

# get the regions of the associated trans genes
trans_genes <- eqtl$GeneSymbol
trans_genes <- gene_annot[gene_annot$SYMBOL %in% trans_genes]

# ------------------------------------------------------------------------------
print("Retrieving SNP genes.")
# ------------------------------------------------------------------------------

# get the relevant snp genes by overlapping with our sentinel region
snp_genes <- subsetByOverlaps(promoters(gene_annot),
                               sentinel_extrange,
                               ignore.strand=T)

# ------------------------------------------------------------------------------
print("Collecting TFs and shortest path genes.")
# ------------------------------------------------------------------------------

# load all TFBS we have available in our data and connect with trans-genes
tfbs <- tfbs[trans_genes$SYMBOL,,drop=F]

# init
tfs <- NULL
sp <- NULL

# get TFs and map their corresponding trans gene
tfs_by_transGene <- c()
for(i in 1:length(trans_genes)){
  s <- trans_genes[i]$SYMBOL
  tfbs_sub <- tfbs[s,,drop=F]
  tfbs_sub <- names(tfbs_sub[,apply(tfbs_sub,2,any)])
  if(length(tfbs_sub)>0) {
    tfbs_sub <- unique(unlist(lapply(strsplit(tfbs_sub, "\\."), "[[", 1)))
    tfs <- gene_annot[gene_annot$SYMBOL %in% tfbs_sub]
    n <- names(tfs_by_transGene)
    tfs_by_transGene <- c(tfs_by_transGene, tfs)
    names(tfs_by_transGene) <- c(n,s)
  }
}
tfs <- unique(unlist(GenomicRangesList(tfs_by_transGene)))

# find the shortest path genes between the SNP genes and the annotated TFs
if(length(tfs)<1){
  warning("No TFs, skipping shortest paths calculation.")
} else {
  # add any TFs, trans genes and their binding information to the PPI db
  toadd <- setdiff(tfs$SYMBOL, ppi_genes)
  toadd <- c(toadd, setdiff(trans_genes$SYMBOL, ppi_genes))
  ppi_db_mod <- addNode(toadd, ppi_db)
  for(i in 1:length(tfs_by_transGene)) {
    tf_sub <- tfs_by_transGene[[i]]
    tgene <- names(tfs_by_transGene)[i]
    ppi_db_mod <- addEdge(rep(tgene, length(tf_sub)),
                          tf_sub$SYMBOL,
                          ppi_db_mod)
  }
  syms_sp <- get_shortest_paths(cis = unique(c(tfs$SYMBOL, snp_genes$SYMBOL)),
                                trans=trans_genes$SYMBOL,
                                snp_genes$SYMBOL,
                                ppi_db_mod)

  # did we find any?
  if(length(syms_sp) < 1){
    warning("No shortest path genes.")
  } else {
    sp <- gene_annot[gene_annot$SYMBOL %in% syms_sp]
  }
}

# ------------------------------------------------------------------------------
print("Finalizing and saving results.")
# ------------------------------------------------------------------------------
result <-  list(sentinel = sentinel_range,
                snp_genes=snp_genes,
                trans_genes=trans_genes)

if(!is.null(sp)){
  result$spath <- sp
}
if(!is.null(tfs)){
  result$tfs <- tfs
  result$tfs_by_transGene <- tfs_by_transGene
}

saveRDS(file=fout, result)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
