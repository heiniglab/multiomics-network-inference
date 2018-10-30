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
mediation <- function(data, snp, genes, cpgs, fout=NULL) {
  require(ggplot2)

  if(any(!c(snp,genes,cpgs) %in% colnames(data))){
    stop("Some of the provided entities are not contained in data.")
  }
  if(any(!grepl("^cg", cpgs))){
    warning("CpG ids don't look like probe ids, continuing anyway.")
  }

  d <- data[,c(snp,genes,cpgs),drop=F]

  # get large matrix of all coefficients for all combinations
  betas <- c()
  for(g in genes) {
    if(grepl("-", g)){
      g <- paste0("`",g,"`")
    }

    # snp-gene
    snp.gene <- lm(paste0(g, "~",snp), data=d)
    snp.gene <- coefficients(snp.gene)[snp]

    # snp-cpg and gene-cpg
    for(c in cpgs) {

      snp.cpg <- lm(paste0(c,"~",snp), data=d)
      snp.cpg <- coefficients(snp.cpg)[snp]

      gene.cpg <- lm(paste0(c, "~", g), data=d)
      gene.cpg <- coefficients(gene.cpg)[g]

      # snp-cpg hat
      snp.cpg.hat <- snp.gene * gene.cpg

      # save result
      betas <- rbind(betas, c(snp, g, c,
			      snp.cpg, snp.gene, gene.cpg, snp.cpg.hat))
    }
  }
  colnames(betas) <- c("snp", "gene", "cpg", "snp.cpg", "snp.gene", "gene.cpg", "snp.cpg.hat")
  betas <- as.data.frame(betas, stringsAsFactors=F)
  betas$snp.cpg <- as.numeric(betas$snp.cpg)
  betas$snp.gene <- as.numeric(betas$snp.gene)
  betas$gene.cpg <- as.numeric(betas$gene.cpg)
  betas$snp.cpg.hat <- as.numeric(betas$snp.cpg.hat)

  # plot the betas for each gene
  if(!is.null(fout)) {
    gp <- ggplot(data=betas, aes(x=snp.cpg, y=snp.cpg.hat)) +
	    geom_hline(yintercept=0, colour="grey") +
	    geom_vline(xintercept=0, colour="grey") +
	    facet_wrap(~ gene, ncol=3) +
  	    geom_point() +
  	    geom_smooth(method="lm") +
  	    ylab(expression(hat(beta)[c])) +
  	    xlab(expression(beta[c])) +
  	    geom_abline(slope=1, colour="orange")
    ggsave(fout, plot=gp)
  }

  # get the correlation results for the estimated against the
  # observed betas
  result <- lapply(genes, function(g){
    if(grepl("-", g)){
      g <- paste0("`",g,"`")
    }

    d <- betas[betas$gene == g,]

    d1 <- d[,"snp.cpg"]
    d2 <- d[,"snp.cpg.hat"]

    fit <- lm(d2~d1)
    r <- cor.test(d1, d2)
    cor <- r$estimate
    pv <- r$p.value

    # also construct a contingency table where we only look at the signs
    # of the individual beta coefficients for each CpG and test the association
    tab <- matrix(0,ncol=2,nrow=2)
    for(i in 1:nrow(d)) {
       obs <- d[i,"snp.cpg"]
       pred <- d[i,"snp.cpg.hat"]
       if(obs<0&pred<0) tab[2,2] <- tab[2,2] + 1
       if(obs<0&pred>0) tab[1,2] <- tab[1,2] + 1
       if(obs>0&pred<0) tab[2,1] <- tab[2,1] + 1
       if(obs>0&pred>0) tab[1,1] <- tab[1,1] + 1
    }
    ft <- fisher.test(tab)
    pv_fisher <- ft$p.value
    or <- ft$estimate

    c(correlation=unname(cor),
      pval_cor=pv,
      odds_ratio=unname(or),
      pval_fisher=pv_fisher)
  })
  names(result) <- genes
  return(result)
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
#' @param cutoff Significance cutoff for mediating genes
#'
#' @author Johann Hawe
#'
mediation.summary <- function(med, s, s.selected, cutoff) {
  # mediation for only ggm selected snp genes
  med.selected <- med[s.selected]
  med.selected.sign <- med.selected[unlist(lapply(med.selected,
                                                  function(x) {
                                                    x["pval_fisher"]<cutoff
                                                  }))]

  # mediation for only not selected snp genes
  med.notselected <- med[setdiff(s,s.selected)]
  med.notselected.sign <- med.notselected[unlist(lapply(med.notselected,
                                function(x) {
                                  x["pval_fisher"]<cutoff
                                }))]

  # get the best mediating gene based on the correlation
  best <- med[which.max(abs(unlist(lapply(med, "[[", "odds_ratio"))))]
  best_mediator <- names(best)
  best_odds_ratio <- unname(best[[1]]["odds_ratio"])

  # TODO if we choose to do a test for comparing
  # selected/not selected, check assumptions first

  # min mediation pv for all not selected snpgenes
  med.pv <- min(unlist(lapply(med.notselected, "[[", "pval_fisher")))
  # max mediation pv for selected snpgenes
  med.pv.selected <- max(unlist(lapply(med.selected, "[[", "pval_fisher")))

  cat("Mediation result summary (min/max of pvals):\n")
  cat("not selected\tselected\n")
  cat(med.pv, "\t", med.pv.selected, "\n")

  # compare the two pvalues
  d <- log10(med.pv/med.pv.selected)
  cat("Difference (log10(ns/s)): ", d, "\n")

  # return 'nice' result
  return(c(length(which(unlist(med)<cutoff)),
	   best_mediator,
	   best_odds_ratio,
	   length(med.notselected.sign),
           length(med.selected.sign),
           paste(names(med.notselected.sign), collapse=";"),
           paste(names(med.selected.sign), collapse=";"),
           paste(names(med.notselected), collapse=";"),
           paste(names(med.selected), collapse=";"),
           paste(format(med.notselected, digits=3), collapse=";"),
           paste(format(med.selected, digits=3), collapse=";"),
           med.pv,
           med.pv.selected,
           d))
}

#' Creates a summary of the comparison of mediation results
#' from two different datasets
#'
#' Given mediation results for a SNP and it's SNP-genes together
#' with the GGM selected SNP genes, retrieves some values indicating
#' the mediation 'performance' for the current SNP
#'
#' @param sentinel The corresponding sentinel snp
#' @param mediation The mediation results on the validation cohort
#' @param mediation2 The mediation results on the original cohort
#' @param sgenes_selected List of all snp genes selected by the model
#' @param cutoff Significance cutoff for mediating genes
#' @param fout Output file were plot should be saved to
#'
#' @author Johann Hawe
#'
compare_mediation_results <- function(sentinel,
				      mediation,
				      mediation2,
				      sgenes_selected,
				      cutoff,
				      fout) {
  require(ggplot2)
  require(reshape)

  # get mediation correlations for all genes
  med_correlations <- unlist(lapply(mediation, "[[", "correlation"))
  med_correlations2 <- unlist(lapply(mediation2, "[[", "correlation"))
  med_correlations2 <- med_correlations2[names(med_correlations)]

  # get mediation pvalues for all genes
  med_pvals <- unlist(lapply(mediation, "[[", "pval_fisher"))
  med_pvals2 <- unlist(lapply(mediation2, "[[", "pval_fisher"))
  med_pvals2 <- med_pvals2[names(med_pvals)]

  # create dataframe for plotting
  df <- cbind.data.frame(med_correlations, med_correlations2,
                          med_pvals, med_pvals2)
  colnames(df) <- c("correlation_validation",
		    "correlation_original",
		    "pvalue_validation",
		    "pvalue_original")

  df$gene <- names(med_pvals)
  df$significant_validation <- df$pvalue_validation <= cutoff
  df$significant_original <- df$pvalue_original <= cutoff
  df$significant <- unlist(lapply(1:nrow(df), function(i) {
			  r <- df[i,,drop=T]
			  if(r[["significant_validation"]] & r[["significant_original"]]){
			    return("both")
			  } else if(r[["significant_original"]]){
			    return("original")
			  } else if(r[["significant_validation"]]) {
			    return("validation")
			  } else {
			    return("none")
			  }
  }))
  df$model_selected <- ifelse(df$gene %in% sgenes_selected, "yes", "no")
  df <- reshape::melt(df, measure.vars=1:2, variable.name="dataset")

  # add plotting x-lab proxy
  df = transform(df, dvariable = as.numeric(variable) + runif(length(variable),-0.2,0.2))

  # create the plot
  source("scripts/lib.R")
  cols <- set_defaultcolors()
  sfm <- scale_fill_manual(values=cols)

  gp <- ggplot(data=df, aes(y=value, x=variable))
  gp <- gp + sfm + geom_violin(draw_quantiles=.5) +
	  geom_point(aes(color=significant, x=dvariable, shape=model_selected), size=3) +
	  geom_line(aes(group=gene, x=dvariable), linetype="dotted")+
	  geom_text(aes(label=gene), size=2.5, nudge_x=-0.4,
		    data=subset(df, variable=="correlation_validation")) +
	  ggtitle(paste0("Mediation correlation values for ", sentinel, "."),
			 "Shown are the distribution of correlation values for each SNP gene, once for the
validation dataset and once for the original dataset (on which models were calcualted.") +
	  theme(plot.title = element_text(hjust = 0))
  ggsave(file=fout,
	       plot=gp)

  # return the correlation of the mediation correlations between the two datasets
  # and the fraction of significant mediations in validation set also significant in original
  # dataset
  corr <- cor(med_correlations, med_correlations2)
  df <- subset(df, variable=="correlation_validation")
  # fraction of all SNP genes significant in both datasets
  frac <- length(which(df$significant=="both"))/nrow(df)
  # fraction of validation significant SNP genes significant in original dataset
  frac2 <- length(which(df$significant=="both")) /
	  length(which(df$significant=="both" | df$significant=="validation"))
  return(c(corr, frac, frac2))
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

  print("Confusion matrix for cpg-genes:")
  cont <- matrix(0, nrow=2, ncol=2)
  cont[1,1] <- v3
  cont[1,2] <- v1
  cont[2,1] <- v4
  cont[2,2] <- v2

  rownames(cont) <- c("has_eQTM", "has_no_eQTM")
  colnames(cont) <- c("ggm_selected", "not_gmm_selected")
  print(cont)

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
    dset <- dset[,colnames(dset) %in% gnodes]
    if(ncol(dset) == 0) {
      return(NA)
    }
    # create adj_matrix
    m <- matrix(nrow=ncol(dset), ncol=ncol(dset))
    rownames(m) <- colnames(m) <- colnames(dset)

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
  				  lambda=seq(0.05,max(pvs$pval), 0.05))$qvalues)
      use <- pvs$qval<0.05 & abs(pvs$cor)>0.3
    } else {
      pvs <- cbind(pvs, qval=rep(NA, nrow(pvs)))
      use <- rep(F, nrow(pvs))
    }

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
validate.geneenrichment <- function(gnodes) {

  # get only gene nodes
  gn <- gnodes[!grepl("^rs|^cg",gnodes)]

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
