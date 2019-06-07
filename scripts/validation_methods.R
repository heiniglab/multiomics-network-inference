#' Library script providing the main functions for validating a ggm.fit based
#' on the identified links between entities
#'
#' @author Johann Hawe
#'
#' @version 20170803
#'


#' Validates the given ggm fit based on the found SNP~Gene links
#'
#' @param graph The graph extracted from the ggm fit
#' @param data The original data with which the ggm was fitted
#' @param ranges The original list of granges related to the entities in the
#' fitted data matrix
#'
#' @return
#'
#' @author Johann Hawe
#'
validate.snps <- function(graph, data, ranges){

}

#' Validates the given ggm fit based on the found CpG~Gene links
#'
#' @param ggm.fit The ggm.fit which is to be validated
#' @param ranges The original list of granges related to the entities in the
#' fitted data matrix
#'
#' @return
#'
#' @author Johann Hawe
#'
validate.cpgs <- function(ggm.fit, ranges){
  # (1) Epigenetic annotation chromHMM
  # (2) trans-eQTL in other dataset (Cpg~Gene)
  # (3) CpG in TFBS in independent CLs
}

#' Validates cpg genes based on a set of eqtm-genes
#'
#' Rather naive approach for validating our identified cpg genes.
#' Checks which of the selected and not selected cpg genes were identified
#' as having an eQTM (given by a list of eqtm.genes) and creates from this
#' a confusion table. In the end the fisher.test pvalue of this table is
#' reported
#'
#' @param eqtm.genes A list of eqtm genes which should be checked against
#' @param cg The list of cpg-genes in our network
#' @param cg.selected The list of cpg-genes in out network selected by the
#' GGM
#'
#' @author Johann Hawe
#'
#' @return The fisher.test()-pvalue for the calculated confusion table
#'
validate_cpggenes <- function(eqtm.genes, cg, cg.selected){

  # create the individual values for our confusion matrix
  v1 <- sum(setdiff(cg, cg.selected) %in% eqtm.genes)
  v2 <- sum(!setdiff(cg, cg.selected) %in% eqtm.genes)
  v3 <- sum(cg.selected %in% eqtm.genes)
  v4 <- sum(!cg.selected %in% eqtm.genes)

  # build the confusion matrix

  print("Confusion matrix for cpg-genes:")
  cont <- matrix(0, nrow=2, ncol=2)
  cont[1,1] <- v3
  cont[1,2] <- v1
  cont[2,1] <- v4
  cont[2,2] <- v2

  rownames(cont) <- c("has_eQTM", "has_no_eQTM")
  colnames(cont) <- c("ggm_selected", "not_gmm_selected")

  # report fisher test
  return(c(bonder_cis_eQTM=fisher.test(cont)$p.value))
}

#' Validates the given ggm fit based on the found Gene~Gene links
#'
#' @param expr.data A list containing expression datasets (matrices) which
#' should be checked against
#' @param g The GGM graph for which to check the genes
#' @param all.genes The list of genes intitially gone into the GGM analysis
#'
#' @author Johann Hawe
#'
#' @return For each of the expression datasets a single (fisher) pvalues for the
#' respective set
#'
validate_gene2gene <- function(expr.data, g, all.genes){

  require(BDgraph)

  # the complete set of nodes in the ggm graph
  gnodes <- nodes(g)

  # get adjacency matrices
  model_adj <- as(g, "matrix")

  # we create a adjacency matrix for each of the expression sets
  # TODO: for the external datasets, check whether we could normalize for age/sex
  results <- lapply(names(expr.data), function(ds) {
    # get the data set
    dset <- expr.data[[ds]]
    dset <- dset[,colnames(dset) %in% gnodes,drop=F]
    if(ncol(dset) < 2) {
      return(NA)
    }
    # create adj_matrix
    m <- matrix(nrow=ncol(dset), ncol=ncol(dset))
    rownames(m) <- colnames(m) <- colnames(dset)
    diag(m) <- 0

    # calculate correlations
    res <- c()
    cnames <- colnames(dset)
    for(i in 1:ncol(dset)) {
      for(j in i:ncol(dset)) {
      	if(j==i) next
        corr <- cor.test(dset[,i], dset[,j])
        res <- rbind(res, c(cnames[i], cnames[j], corr$p.value, corr$estimate))
      }
    }
    pvs <- as.data.frame(matrix(res, ncol=4, byrow = F), stringsAsFactors=F)
    colnames(pvs) <- c("n1", "n2", "pval", "cor")
    pvs$pval <- as.numeric(pvs$pval)
    pvs$cor <- as.numeric(pvs$cor)

    # get qvalue
    if(nrow(pvs)>10) {

      pvs <- cbind(pvs, qval=qvalue(pvs$pval,
  				  lambda=seq(min(pvs$pval), max(pvs$pval)))$qvalues)
    } else {
      pvs <- cbind(pvs, qval=p.adjust(pvs, "BH"))
    }
    use <- pvs$qval<0.05 & abs(pvs$cor)>0.3
    # fill matrix
    for(i in 1:nrow(pvs)) {
     r <- pvs[i,,drop=F]
     if(use[i]) {
       m[r$n1, r$n2] <- m[r$n2, r$n1] <- 1
     } else {
       m[r$n1, r$n2] <- m[r$n2, r$n1] <- 0
     }
    }
    n <- intersect(colnames(m), colnames(model_adj))
    m <- m[n,n]
    model_adj <- model_adj[n,n]
    res <- BDgraph::compare(model_adj, m)
    res["MCC","estimate"]
  })
  names(results) <- names(expr.data)
  # report MCC for each dataset
  return(unlist(results))
}

#' Validate the genes in the GGM by performing a
#' GO enrichment on them
#'
#' @param gnodes The nodes in the GGM graph
#'
#' @author Johann Hawe
#'
#' @return A vector containing enriched GOIDs, terms, their pvalues and
#' qvalues or instead only containing NAs if no enrichments was found
#'
validate_geneenrichment <- function(gnodes) {

  # get only gene nodes
  gn <- gnodes[!grepl("^rs|^cg",gnodes)]

  if(length(gn) < 2) return(NULL)

  # define background set
  # for now all possible symbols from the array annotation
  # are used
  #bgset <- dnodes[!grepl("^rs|^cg",dnodes)]

  annot <- illuminaHumanv3SYMBOLREANNOTATED
  bgset <- unique(unlist(as.list(annot)))

  if(length(gn)>0 & length(gn)<length(bgset)){
    go.tab <- go.enrichment(gn, bgset, gsc)
    go.tab <- go.tab[go.tab$q<0.01,,drop=F]
    if(nrow(go.tab)>0){
      r <- c(paste0(go.tab$GOID, collapse=","),
               paste0(go.tab$Term, collapse=","),
               paste0(go.tab$Pvalue, collapse=","),
               paste0(go.tab$q, collapse=","))
      return(r)
    }
  }
  r <- c(NA, NA, NA, NA)
  return(r)
}

#' Validates all trans genes (tfs and cpg-genes) by using a set of
#' previously identified trans eqtls
#'
#' @param teqtl The set of related trans eqtls
#' @param cgenes The CpG-genes to be checked
#' @param cgenes.selected THe CpG-genes selected in the GGM
#' @param tfs The TFs to be checked
#' @param tfs.selected The TFs selected in the GGM
#'
#' @author Johann Hawe
#'
#' @return Fisher pvalues for the two contingency table tests
#'
validate_trans_genes <- function(teqtl, cgenes, tfs,
                                 cgenes.selected, tfs.selected) {

  # analyze the cpggenes, total and selected
  teqtl.cgenes <- cgenes[cgenes %in% unlist(strsplit(teqtl$Transcript_GeneSymbol, "\\|"))]
  teqtl.cgenes.selected <- intersect(teqtl.cgenes, cgenes.selected)

  # analyze the tfs, total and selected
  teqtl.tfs <- tfs[tfs %in% unlist(strsplit(teqtl$Transcript_GeneSymbol, "\\|"))]
  teqtl.tfs.selected <- intersect(teqtl.tfs, tfs.selected)

  cat("CpG genes: ", cgenes, "\n")
  cat("CpG genes with trans-eQTL: ", teqtl.cgenes, "\n")
  cat("selected CpG genes with trans-eQTL: ", teqtl.cgenes.selected, "\n")

  cat("TFs: ", tfs, "\n")
  cat("TFs with trans-eQTL: ", teqtl.tfs, "\n")
  cat("selected TFs with trans-eQTL: ", teqtl.tfs.selected, "\n")

  # create matrix for fisher test
  cont <- matrix(c(length(teqtl.cgenes),length(teqtl.cgenes.selected),
                   length(cgenes),length(cgenes.selected)),
                 nrow=2,ncol=2, byrow = T)
  cat("confusion matrix for cgenes:\n")
  rownames(cont) <- c("teqtl", "no teqtl")
  colnames(cont) <- c("not selected", "selected")
  f1 <- fisher.test(cont)$p.value

  # create matrix for fisher test
  cont <- matrix(c(length(teqtl.tfs),length(teqtl.tfs.selected),
                   length(tfs),length(tfs.selected)),
                 nrow=2,ncol=2, byrow = T)
  cat("confusion matrix for TFs:\n")
  rownames(cont) <- c("teqtl", "no teqtl")
  colnames(cont) <- c("not selected", "selected")
  f2 <- fisher.test(cont)$p.value

  # report fisher test results
  return(c(transEqtl_cgenes=f1,transEqtl_tfs=f2))
}
