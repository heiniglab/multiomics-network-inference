#'
#' Script to collect all associated ids (expr, meth and geno)
#' for a specific sentinel snp
#'
#' @author Johann Hawe
#' 
#' @date 2017/03/04
#'

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
#' @param cis List of nodes
#' @param trans List of nodes
#' @param string_db Instance of the string database to be used (graphNEL)
#' 
#' @return A vector of gene symbols being on the shortest path between the two lists
#' of genes as found in the validated string network
#'
#' @author Johann Hawe, Matthias Heinig
#' 
get_string_shortest_paths <- function(cis, trans, snp_genes, string_db) {
  
  g <- string_db
  nodes <- nodes(g)
  
  # ----------------------------------------------------------------------------
  # ensure to have only nodes in our giant cluster
  # ----------------------------------------------------------------------------
  cis <- cis[which(cis %in% nodes)]
  trans <- trans[which(trans %in% nodes)]
  snp_genes <- snp_genes[which(snp_genes %in% nodes)]
  if(length(cis) == 0 | length(trans) == 0 | length(snp_genes) == 0) {
    return(NULL)
  }
  
  # ----------------------------------------------------------------------------
  # calculate weights for nodes
  print("Starting propagation.")
  # ----------------------------------------------------------------------------
  prop = propagation(graph2sparseMatrix(g), n.eigs=500, 
                     from=cis, to=trans, sum="both")
  node_weights <- rowSums(prop)
  ## find the shortest path with maximal weight
  ## we have an algorithm that finds minimum node weight paths so we need
  ## to turn the weighting around
  node_weights = max(node_weights) - node_weights + 1
  
  # get the best trans gene as target for SP calculation
  best_trans = snp_genes[which.max(prop[snp_genes,"from"])]
  print(paste0("Best trans genes: ", paste(best_trans, collapse=", ")))
  
  # ----------------------------------------------------------------------------
  # extract shortest paths
  print("Extracting shortest paths to best trans genes.")
  # ----------------------------------------------------------------------------
  sp = min.node.weight.path(g, node_weights, from=cis, to=best_trans)
  nodes <- setdiff(unlist(lapply(sp, "[", "path_detail")), NA)
  nodes <- setdiff(nodes, c(trans, cis))
  
  return(nodes)
}

#' Gets all ranges-data for a specific sentinel SNP
#'
#' Check for each sentinel SNP which genes are near it and its associated 
#' CpGs and gets their ranges. Also saves the ranges of the sentinel 
#' and its other associated SNPs. Filters by number of CpGs per sentinel.
#'
#' @param id The id of sentinel to be processed.
#' @param cosmo The individual meqtl associations
#' @param trans_meQTL The pruned trans-meQTL table
#' 
#' @value A list containing a relevant ranges for the sentinel snp:
#'  cpg ranges (meQTL cpgs for sentinel) -> cpgs
#'  snp ranges -> snp_ranges
#'  sentinel ranges -> sentinel_range and sentinel_ext_range
#'  cpg and snp gene ranges -> cpg_genes and snp_genes
#'  tfs -> TFs with TFBS at CpGs
#'  spath -> genes on shortest path between CpGs and snp genes
#'  
#' @author Johann Hawe
#' 
#' @date 02/07/2017
#'
collect_ranges <- function(id, cosmo, trans_meQTL, 
                           fcpgcontext, sentinel_only=T) {
  
  # ----------------------------------------------------------------------------
  print("collecting SNP and CpG ranges.")
  # ----------------------------------------------------------------------------
  
  pairs = which(trans_meQTL[,"sentinel.snp"] == id)
  
  # get sentinel idx
  sentinel_idx <- which(cosmo$snp==id)[1]
  
  # the large interval range for the sentinel
  chr <- paste0("chr", trans_meQTL[pairs,"chr.snp"][1])
  start <- trans_meQTL[pairs,"interval.start.snp"][1]
  end <- trans_meQTL[pairs,"interval.end.snp"][1]
  
  sentinel_extrange <- GRanges(chr, IRanges(start,end))
  sentinel_range <- with(cosmo[sentinel_idx,], 
                         GRanges(paste0("chr", snp.chr),
                         IRanges(snp.pos, width=1)))
  names(sentinel_range) <- names(sentinel_extrange) <- id
  
  if(!sentinel_only) {
    # get ranges of SNP which the sentinel represents
    sentinel_subset <- trans_meQTL[pairs,]
    snp_ranges <- GRanges(paste0("chr", sentinel_subset[,"chr.snp"]), 
                          IRanges(sentinel_subset[,"pos.snp"],width=1))
    names(snp_ranges) <- sentinel_subset[,"snp"]
    # remove the sentinel from the list of additional SNPs
    snp_ranges <- unique(snp_ranges[!grepl(paste0(id,"$"), names(snp_ranges))])
  }
  
  # get cosmo subset
  idxs <- get.trans.cpgs(id, trans_meQTL, cosmo)
  
  # get related genes, i.e. genes near meQTL loci (snp+cpg)
  
  # extend cpgs
  cosmosub <- cosmo[idxs,]
  
  croi <- with(cosmosub, GRanges(paste0("chr", cpg.chr), 
                                 IRanges(cpg.pos,width=2)))
  names(croi) <- as.character(cosmosub[,"cpg"])
  croi <- unique(croi)
  
  # extended sentinel region
  sroi <- sentinel_extrange
  names(sroi) <- id
  
  # ----------------------------------------------------------------------------
  print("Retrieving SNP and CpG genes.")
  # ----------------------------------------------------------------------------
  
  # for now use the results from the enrichment analysis to get the
  # relevant SNP genes
  genes_sroi <- gene_annot[gene_annot %over% sroi]
  genes_sroi$ids <- probes.from.symbols(genes_sroi$SYMBOL, as.list=T)
  if(is.null(genes_sroi$ids)){
    warning(paste0("No expr probe ids available for snp genes for sentinel ",
                   id))
    return(NULL)
  }
  
  # get genes near our cpg regions
  genes_croi_by_cpg <- get.nearby.ranges(croi, promoters(gene_annot))
  temp <- lapply(names(genes_croi_by_cpg), function(cg) {
    prbs <- probes.from.symbols(genes_croi_by_cpg[[cg]]$SYMBOL, as.list=T)
    if(is.null(prbs)) {
      genes_croi_by_cpg[[cg]]$ids <<- NA
    } else {
      if(any(is.null(prbs))) {
        prbs[[unlist(lapply(prbs,is.null))]] <- NA
      }
      genes_croi_by_cpg[[cg]]$ids <<- prbs
    }
    prbs
  })
  
  # get single list of all cpg genes
  genes_croi <- unlist(GRangesList(unlist(genes_croi_by_cpg)))
  
  if(is.null(unlist(temp))){
    warning(paste0("No expr probe ids available for cpg genes for sentinel ",
                   id))
    return(NULL)
  }
  
  # ----------------------------------------------------------------------------
  print("Collecting TFs and shortest path genes.")
  # ----------------------------------------------------------------------------
  
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
  locus_graph <- add.to.graphs(list(string_db), id, snp_genes, 
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
    
    print("Finding shortest path genes.")
    syms_sp <- get_string_shortest_paths(cis = cpgs_with_tfbs, 
                                         trans=unique(c(snp_genes_in_string,
                                                        tfs)), 
                                         snp_genes=snp_genes_in_string,
                                         locus_graph)
    
    print(paste0("Shortest paths genes: ", paste(syms_sp, collapse=", ")))
    
    if(length(syms_sp) < 1){
      warning("No shortest path genes.")
    } else {
      sp <- gene_annot[gene_annot$SYMBOL %in% syms_sp]
      sp$ids <- probes.from.symbols(sp$SYMBOL, as.list=T)
    }
    tfs <- gene_annot[gene_annot$SYMBOL %in% tfs]
    tfs$ids <- probes.from.symbols(tfs$SYMBOL, as.list=T)
  }
  
  # construct our result list
  result <-  list(cpgs=croi,sentinel_range=sentinel_range,
                  sentinel_ext_range=sentinel_extrange, 
                  snp_genes=genes_sroi, cpg_genes=genes_croi,
                  cpg_genes_by_cpg=genes_croi_by_cpg)
  
  if(!sentinel_only) {
    result$snp_ranges <- snp_ranges
  }
  if(!is.null(sp)){
    result$spath <- sp
  }
  if(!is.null(tfs)){
    result$tfs <- tfs
  }
  return(result)
}

# ------------------------------------------------------------------------------
# Start processing
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------
fcosmo <- snakemake@input[["tcosmo"]]
fmeqtl <- snakemake@input[["meqtl"]]
fstring <- snakemake@input[["string"]]

# TODO: create this file from scratch!
fcpgcontext <- snakemake@input[["cpgcontext"]]

ofile <- snakemake@output[[1]]
sentinel <- snakemake@params[["sentinel"]]

# ------------------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------------------

cat("Loading data.\n")

gene_annot <- get.gene.annotation()
string_db <- readRDS(fstring)
  
# load trans-meQTL table
trans_meQTL = read.csv(fmeqtl, 
                       sep="\t", stringsAsFactors=F)

# load trans-cosmo information
cosmo <- readRDS(fcosmo)

# ------------------------------------------------------------------------------
# Collect and save ranges
# ------------------------------------------------------------------------------

print("Collecting ranges.")

# collect all associated ranges
ranges <- collect.ranges(sentinel, cosmo, trans_meQTL, fcpgcontext)

saveRDS(file=ofile, ranges)
