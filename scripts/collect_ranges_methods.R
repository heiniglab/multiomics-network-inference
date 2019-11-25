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
#' @param cis List of nodes in cis (e.g. CpGs)
#' @param trans List of nodes in trans (e.g. SNP genes)
#' @param snp_genes List of snp genes
#' @param ppi_db Instance of the PPI database to be used (graphNEL)
#' @param tfs The TFs (symbols) collected for the locus
#' @param best_trans Optional: previously identified 'best' trans gene (e.g.
#' according to a random walk as in the meQTL paper.
#'
#' @return A vector of gene symbols being on the shortest path between the
#' two lists of genes as found in the validated string network
#'
#' @author Johann Hawe, Matthias Heinig
#'
# ------------------------------------------------------------------------------
get_shortest_paths <- function(cis, trans, snp_genes, ppi_db, 
                               tfs, best_trans=NULL) {

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
  sp = min.node.weight.path(ppi_db, node.weights, from=cis, to=best_snp_gene)
  nodes <- setdiff(unlist(lapply(sp, "[", "path_detail")), NA)
  nodes_non_tf <- setdiff(nodes, c(trans, cis))
  nodes_tf <- intersect(nodes, tfs)

  print("Shortest path genes:")
  print(nodes_non_tf)
  print("TFs on shortest paths:")
  print(nodes_tf)

  return(list(non_tf_sp=nodes_non_tf, 
              tf_sp=nodes_tf))
}

#' ------------------------------------------------------------------------------
#' Get for trans genes the TFs which bind in their promoter
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' ------------------------------------------------------------------------------
get_tfs_by_transGene <- function(tfbs, trans_genes, gene_annot) {
  # get TFs and map their corresponding trans gene
  trans_genes <- trans_genes[trans_genes$SYMBOL %in% rownames(tfbs)]
  print("Trans genes:")
  print(trans_genes)
  tfs_by_transGene <- c()
  for(i in 1:length(trans_genes)){
    s <- trans_genes[i]$SYMBOL
    n <- names(trans_genes)[i]
    tfbs_sub <- tfbs[s,,drop=F]
    tfbs_sub <- unique(colnames(tfbs_sub[,apply(tfbs_sub,2,function(x) sum(x) > 0),drop=F]))
    if(length(tfbs_sub)>0) {
      tf_symbols <- unique(sapply(strsplit(tfbs_sub, "\\."), "[[", 1))
      tfs <- gene_annot[gene_annot$SYMBOL %in% tf_symbols]
      na <- names(tfs_by_transGene)
      tfs_by_transGene <- c(tfs_by_transGene, unique(tfs))
      names(tfs_by_transGene) <- c(na,s)
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
                                        ppi_genes, snp_genes, ppi_db,
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
  shortest_paths <- get_shortest_paths(cis = trans_genes,
                                trans=unique(c(snp_genes, tfs)),
                                snp_genes,
                                ppi_db_mod,
                                tfs)

  non_tf_sp <- shortest_paths$non_tf_sp
  tf_sp <- shortest_paths$tf_sp

  # did we find any?
  if(length(non_tf_sp) < 1){
    warning("No shortest path genes.")
    non_tf_sp_ranges <- NULL
  } else {
    non_tf_sp_ranges <- gene_annot[gene_annot$SYMBOL %in% non_tf_sp]
  }
  # check also TFs on shortest paths
  if(length(tf_sp) < 1){
    stop("No TF on the shortest paths.")
  } else {
    tf_sp_ranges <- gene_annot[gene_annot$SYMBOL %in% tf_sp]
  }

  return(list(non_tf_sp = non_tf_sp_ranges,
              tf_sp = tf_sp_ranges))
}

# ------------------------------------------------------------------------------
# Simple helper to quickly plot ranges
# ------------------------------------------------------------------------------
plot_ranges <- function(ranges, fout) {
  require(ggplot2)
  cols <- get_defaultcolors()

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

  pdf(fout, width=12, height=12)

  # barplot of number of individual entity types
  print(ggplot(aes(y=count, x=name), data=toplot) +
    geom_bar(stat="identity", position="dodge") + gt + theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)))

  dev.off()
}

#' -----------------------------------------------------------------------------
#' Gets the first ranges in subject, which are up-/down-stream and overlapping
#' a range in the query
#'
#' @param query ranges for which to get nearby genes
#' @param subject ranges in which to look for nearby genes
#'
#' @return Either the idx of the hits in the subject if idxs=T, or the identified ranges (idxs=F)
#'
#' -----------------------------------------------------------------------------
get.nearby.ranges <- function(query, subject) {

  nearby <- function(q,s){
    # get preceding, following and overlapping instances of any
    # range in query within subject ranges
    pre <- precede(q, s, select="all", ignore.strand=T)
    fol <- follow(q, s, select="all", ignore.strand=T)
    ove <- findOverlaps(q, s, select="all", ignore.strand=T)
    # combine hits
    h <- unique(c(subjectHits(pre), subjectHits(fol), subjectHits(ove)))

    return(list(hits=h, ranges=s[h]))
  }

  #return a list, where each query gets its nearby ranges annotated with their distance
  res <- lapply(1:length(query), function(j) {
    q <- query[j]
    n <- nearby(q, subject)
    if(length(n$hits)>0){
      n$ranges$distance <- rep(-1, times=length(n$ranges))
      n$ranges$hit_idx <- rep(NA_integer_, times=length(n$ranges))
      for(i in 1:length(n$ranges)) {
        d <- distance(q,n$ranges[i])
        n$ranges[i]$distance <- d
        n$ranges[i]$hit_idx <- n$hits[i]
      }
      return(n$ranges)
    } else {
      return(NULL)
    }
  })
  return(res)
}
