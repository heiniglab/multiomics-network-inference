#' Library script providing the main functions for validating a ggm.fit based
#' on the identified links between entities
#' 
#' @author Johann Hawe
#'
#' @version 20170803
#' 


#' Performs mediation analysis for a given set of SNPs, Genes
#' and CpGs.
#'
#' @param data data matrix with individuals in the rows and
#' SNPs, CpGs and Gene entries in the columns
#'
#' @param snps SNP ids to be analyzed
#' @param genes The genes to be checked
#' @param cpgs the cpgs to be checked
#' @param plot Flag whether to create a correlation plot for each of the genes
#' 
#' @author Johann Hawe
#' 
#' 
mediation <- function(data, snp, genes, cpgs, plot=F) {
  if(any(!c(snp,genes,cpgs) %in% colnames(data))){
    stop("Some of the provided entities are not contained in data.")
  }
  if(any(!grepl("^cg", cpgs))){
    warning("CpG ids don't look like probe ids, continuing anyway.")
  }
  
  d <- data[,c(snp,genes,cpgs),drop=F]
  
  # calculate all linear models and get the betas
  # snp-gene
  snp.gene <- lapply(genes, function(g){
    r <- lm(paste0("`", g, "`~",snp), data=d)
    coefficients(r)[snp]
  })
  names(snp.gene) <- genes
  
  # snp-cpg
  snp.cpg <- lapply(cpgs, function(c){
    r <- lm(paste0(c,"~",snp), data=d)
    coefficients(r)[snp]
  })
  names(snp.cpg) <- cpgs
  
  betas <- lapply(genes, function(g){
    # gene-cpg
    temp <- lapply(cpgs, function(c){
      r <- lm(paste0(c, "~`", g, "`"), data=d)
      # merge with the other betas (for each cpg and gene combi, we have 3 betas)
      r <- c(snp.cpg[[c]], 
             snp.gene[[g]], 
             coefficients(r)[2])
      # the estimated beta for snp-cpg via gene
      r <- c(r, r[2] * r[3])
      names(r) <- c(c,g,paste0(c,".",g), paste0(c, ".", g, ".hat"))
      r
    })
  })
  
  betas <- matrix(unlist(betas), 
                  ncol=4,
                  nrow = length(unlist(betas))/4, 
                  byrow = T)
  rnames <- rep("", nrow(betas))
  for(j in 1:length(cpgs)) {
    for(i in 1:length(genes)) {
      n <- paste(genes[i], cpgs[j], sep="."); 
      rnames[i*j+(i-1)*(length(cpgs)-j)] <- n
    }
  }
  rownames(betas) <- rnames
  colnames(betas) <- c("snp.cpg", "snp.gene", "cpg.gene", "snp.cpg.hat")
  
  betas <- as.data.frame(betas)
  betas$cpg <- gsub(".*[.]", "", rownames(betas))
  betas$gene <- gsub("[.].*", "", rownames(betas))
  result <- lapply(genes, function(g){
    d <- betas[betas$gene == g,]
    
    d1 <- d[,"snp.cpg"]
    d2 <- d[,"snp.cpg.hat"]
    fit <- lm(d2~d1)
    r <- cor.test(d1, d2)
    cor <- r$estimate
    pv <- r$p.value
    
    if(plot){
      pdf(paste0("results/current/plots/", snp, "_", g, "_mediation.pdf"))
      plot(x = d1,y = d2,
           main=paste0(g),
           ylab=expression(hat(beta)[c]),
           xlab=expression(beta[c]),
           pch="+")
      abline(fit, col="red")
      text(x=-0.01, y=max(d1), labels = paste0("cor=", format(cor, digits=2), 
                             "\np.value=", format(pv, digits=2)))
      dev.off()
    }
    c(correlation=unname(cor),pvalue=pv)
  })
  names(result) <- genes
  result
}

#' Creates a summary of the mediation results
#' 
#' Given mediation results for a SNP and it's SNP-genes together
#' with the GGM selected SNP genes, retrieves some values indicating
#' the mediation 'performance' for the current SNP
#' 
#' @param m The mediation result gotten from mediation()
#' @param s The SNP-genes for the current instance
#' @param s.selected The SNP-genes selected by the GGM
#' 
#' @author Johann Hawe
#' 
mediation.summary <- function(med, s, s.selected, cutoff=0.05) {
  # mediation for only ggm selected snp genes
  med.selected <- med[s.selected]
  med.selected.sign <- med.selected[unlist(lapply(med.selected, 
                                                  function(x) { 
                                                    x["pvalue"]<cutoff
                                                  }))]
  
  # mediation for only not selected snp genes
  med <- med[setdiff(s,s.selected)]
  med.sign <- med[unlist(lapply(med, 
                                function(x) { 
                                  x["pvalue"]<cutoff
                                }))]
  
  # TODO if we choose to do a test for comparing 
  # selected/not selected, check assumptions first
  
  # min mediation pv for all not selected snpgenes
  med.pv <- min(unlist(lapply(med, "[[", "pvalue")))
  # max mediation pv for selected snpgenes
  med.pv.selected <- max(unlist(lapply(med.selected, "[[", "pvalue")))
  
  cat("Mediation result summary (min/max of pvals):\n")
  cat("not selected\tselected\n")
  cat(med.pv, "\t", med.pv.selected, "\n")
  
  # compare the two pvalues
  d <- log10(med.pv/med.pv.selected)
  cat("Difference (log10(ns/s)): ", d, "\n")
  
  # return 'nice' result
  return(c(length(which(unlist(med)<cutoff)),
	   length(med.sign),
           length(med.selected.sign),
           paste(names(med.sign), collapse=";"),
           paste(names(med.selected.sign), collapse=";"),
           med.pv, 
           med.pv.selected,
           d))
}

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
validate.cpggenes <- function(eqtm.genes, cg, cg.selected){
  
  # create the individual values for our confusion matrix
  v1 <- sum(setdiff(cg, cg.selected) %in% eqtm.genes)
  v2 <- sum(!setdiff(cg, cg.selected) %in% eqtm.genes)
  v3 <- sum(cg.selected %in% eqtm.genes)
  v4 <- sum(!cg.selected %in% eqtm.genes)
  
  # build the confusion matrix 
  cat("confusion matrix for cpg-genes:\n")
  cont <- matrix(c(v3,v1,v4,v2), 
                 ncol=2, byrow=T)
  rownames(cont) <- c("has_eQTM", "has_no_eQTM")
  colnames(cont) <- c("ggm_selected", "not_gmm_selected")
  
  # report fisher test
  return(fisher.test(cont)$p.value)
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
validate.gene2gene <- function(expr.data, g, all.genes){
  
  # the complete set of nodes in the ggm graph
  gnodes <- nodes(g)
  
  # get the gene-gene links which are possible in our graph
  em <- expand.grid(all.genes, all.genes, stringsAsFactors = F)
  em <- (em[em$Var1 != em$Var2, ])
  # remove duplicated entries
  em.names <- vector(mode="character", length=nrow(em))
  for(i in 1:nrow(em)){
    em.names[i] <- paste0(sort(c(em[i,1], em[i,2])), collapse="")
  }
  em <- em[!duplicated(em.names),]
  # create final dataframe to be used
  em <- cbind(em[,1], em[,2])
  em <- as.data.frame(em, stringsAsFactors=F)
  
  # get the gene-gene links which can be found in our graph
  em.selected <- t(graph::edgeMatrix(g))
  em.selected <- cbind(gnodes[em.selected[,1]], gnodes[em.selected[,2]])
  em.selected <- em.selected[!grepl("^rs|^cg", em.selected[,1]),,drop=F]
  em.selected <- em.selected[!grepl("^rs|^cg", em.selected[,2]),,drop=F]
  colnames(em) <- colnames(em.selected) <- c("n1", "n2")
  em.selected <- as.data.frame(em.selected, stringsAsFactors=F)
  
  # we collect for each data set the results (fisher pvalue)
  # TODO: for the external datasets, check whether we could normalize for age/sex
  results <- lapply(names(expr.data), function(ds) {
    # get the data set
    dset <- expr.data[[ds]]
    
    # calculate correlations
    pvs <- lapply(1:nrow(em), function(r){
      n1 <- em[r,1]
      n2 <- em[r,2]
      # check whether we have data for both nodes
      if(n1 %in% colnames(dset) & n2 %in% colnames(dset)){
        res <- cor.test(dset[,n1], dset[,n2])
        c(n1,n2,res$p.value, res$estimate)
      } else {
        c(n1,n2,NA,NA)
      }
    })
    pvs <- as.data.frame(matrix(unlist(pvs), ncol=4, byrow = T), stringsAsFactors=F)
    colnames(pvs) <- c("n1", "n2", "pval", "cor")
    pvs$pval <- as.numeric(pvs$pval)
    pvs$cor <- as.numeric(pvs$cor)
    # get qvalue
    pvs <- cbind(pvs, qval=qvalue(pvs$pval)$qvalues)
    sign.cors <- pvs[pvs$qval<0.01 & abs(pvs$cor)>0.3,]
    
    em.sign <- em[(em$n1 %in% sign.cors$n1 & em$n2 %in% sign.cors$n2)
                  | (em$n2 %in% sign.cors$n1 & em$n1 %in% sign.cors$n2),]
    em.selected.sign <- em.selected[(em.selected$n1 %in% sign.cors$n1 & em.selected$n2 %in% sign.cors$n2)
                                    | (em.selected$n2 %in% sign.cors$n1 & em.selected$n1 %in% sign.cors$n2),]
    
    # create matrix for fisher test
    cont <- matrix(c(nrow(em.selected.sign),nrow(em.sign),
                     nrow(em.selected),nrow(em)),
                   nrow=2,ncol=2, byrow = T)
    rownames(cont) <- c("sign", "not sign")
    colnames(cont) <- c("selected", "not selected")
    print(ds)
    print(cont)
    format(fisher.test(cont)$p.value, digits=4)
  })
  names(results) <- names(expr.data)
  # report fisher test result for each DS
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
validate.geneenrichment <- function(gnodes) {
  
  # get only gene nodes
  gn <- gnodes[!grepl("^rs|^cg",gnodes)]
  
  # define background set
  # for now all possible symbols from the array annotation
  # are used
  #bgset <- dnodes[!grepl("^rs|^cg",dnodes)]
  
  annot <- illuminaHumanv3SYMBOLREANNOTATED
  bgset <- unique(unlist(as.list(annot)))
  
  if(length(gn)<length(bgset)){
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
validate.trans.genes <- function(teqtl, cgenes, tfs,
                                 cgenes.selected, tfs.selected) {
  
  # analyze the cpggenes, total and selected
  teqtl.cgenes <- cgenes[cgenes %in% teqtl$Transcript_GeneSymbol]
  teqtl.cgenes.selected <- intersect(teqtl.cgenes, cgenes.selected)
  
  # analyze the tfs, total and selected
  teqtl.tfs <- tfs[tfs %in% teqtl$Transcript_GeneSymbol]
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
  cont
  f1 <- fisher.test(cont)$p.value
  
  # create matrix for fisher test
  cont <- matrix(c(length(teqtl.tfs),length(teqtl.tfs.selected),
                   length(tfs),length(tfs.selected)),
                 nrow=2,ncol=2, byrow = T)
  cat("confusion matrix for TFs:\n")
  rownames(cont) <- c("teqtl", "no teqtl")
  colnames(cont) <- c("not selected", "selected")
  cont
  f2 <- fisher.test(cont)$p.value

  # report fisher test results
  return(c(f1,f2))
}

sigma.trace <- function(ggm.fit) {
  
}
