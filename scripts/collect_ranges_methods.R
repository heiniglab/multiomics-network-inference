#' ------------------------------------------------------------------------------
#' Define some methods used in the collect ranges scripts.
#'
#' Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' ------------------------------------------------------------------------------

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
# ------------------------------------------------------------------------------
get_shortest_paths <- function(cis, trans, snp_genes, ppi_db, best_trans=NULL) {

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
  if(is.null(best_trans) || !best_trans %in% ppi_genes) {
    best_snp_gene = snp_genes[which.max(prop[snp_genes,"from"])]
  } else {
    best_snp_gene <- best_trans
  }
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

#' ------------------------------------------------------------------------------
#' Get for trans genes the TFs which bind in their promoter
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' ------------------------------------------------------------------------------
get_tfs_by_transGene <- function(tfbs, trans_genes, gene_annot) {
  # get TFs and map their corresponding trans gene
  trans_genes <- trans_genes[trans_genes$SYMBOL %in% rownames(tfbs)]
  tfs_by_transGene <- c()
  for(i in 1:length(trans_genes)){
    s <- trans_genes[i]$SYMBOL
    tfbs_sub <- tfbs[s,,drop=F]
    tfbs_sub <- unique(colnames(tfbs_sub[,apply(tfbs_sub,2,any),drop=F]))
    if(length(tfbs_sub)>0) {
      tfbs_sub <- unique(unlist(lapply(strsplit(tfbs_sub, "\\."), "[[", 1)))
      tfs <- gene_annot[gene_annot$SYMBOL %in% tfbs_sub]
      n <- names(tfs_by_transGene)
      tfs_by_transGene <- c(tfs_by_transGene, unique(tfs))
      names(tfs_by_transGene) <- c(n,s)
    }
  }
  return(tfs_by_transGene)
}

#' ------------------------------------------------------------------------------
#' Get the genes on the shortest path between the trans and cis genes.
#' For all entities, the gene symbols need to be provided here.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' ------------------------------------------------------------------------------
collect_shortest_path_genes <- function(tfs, trans_genes, tfs_by_transGene,
                                        ppi_genes, cis_genes, ppi_db, 
                                        gene_annot) {
  # add any TFs, trans genes and their binding information to the PPI db
  toadd <- setdiff(tfs, ppi_genes)
  toadd <- unique(c(toadd, setdiff(trans_genes, ppi_genes)))
  ppi_db_mod <- addNode(toadd, ppi_db)
  for(i in 1:length(tfs_by_transGene)) {
    tf_sub <- tfs_by_transGene[[i]]
    tgene <- names(tfs_by_transGene)[i]
    ppi_db_mod <- addEdge(rep(tgene, length(tf_sub)),
                          tf_sub$SYMBOL,
                          ppi_db_mod)
  }
  syms_sp <- get_shortest_paths(cis = unique(c(tfs, cis_genes)),
                                trans=trans_genes,
                                cis_genes,
                                ppi_db_mod)
  # did we find any?
  if(length(syms_sp) < 1){
    warning("No shortest path genes.")
    return(NULL)
  } else {
    sp <- gene_annot[gene_annot$SYMBOL %in% syms_sp]
    return(sp)
  }
}

# ------------------------------------------------------------------------------
# Simply helper to quickly plot ranges
# ------------------------------------------------------------------------------
plot_ranges <- function(ranges, fout) {
  require(ggplot2)
  cols <- set_defaultcolors()
  sfm <- scale_fill_manual(values=cols)
  
  # get number of entities
  sg <- length(ranges$snp_genes)
  tfs <- length(ranges$tfs)
  sp <- length(ranges$spath)
  trans_assoc <- length(ranges$trans_genes)

  # create data-frame for use with ggplot
  toplot <- data.frame(name=c("SNP genes", "TFs", "shortest path",
                       "trans associated"),
                       count=c(sg, tfs, sp, trans_assoc))

  gt <- ggtitle(paste0("Entities for sentinel ", names(ranges$sentinel), 
                       " in tissue " , ranges$tissue, "."))

  pdf(fout, width=8, height=8)
  
  # barplot of number of individual entity types
  ggplot(aes(y=count, x=name), data=toplot) +
    geom_bar(stat="identity", position="dodge") + gt +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  dev.off()
}
