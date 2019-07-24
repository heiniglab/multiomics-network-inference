#' -----------------------------------------------------------------------------
#' Methods used in context of the mediation analysis
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jun  4 10:57:48 2019
#' -----------------------------------------------------------------------------

#' -----------------------------------------------------------------------------
#' Performs mediation analysis for a given set of SNPs, Genes
#' and CpGs.
#'
#' @param data data matrix with individuals in the rows and
#' SNPs, CpGs and Gene entries in the columns
#'
#' @param snps SNP ids to be analyzed
#' @param cis_genes The cis/snp genes to be checked
#' @param trans_assoc the trans associated entities to be checked. Either CpG
#' IDs or gene symbols
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
mediation <- function(data, snp, cis_genes, trans_assoc,
                      fout_table, fout_plot) {
  require(ggplot2)

  if(any(!c(snp,cis_genes,trans_assoc) %in% colnames(data))){
    stop("Some of the provided entities are not in the data.")
  }

  d <- data[,c(snp,cis_genes,trans_assoc),drop=F]

  print(length(cis_genes))
  # get large matrix of all coefficients for all combinations
  betas <- c()
  for(g in cis_genes) {
    if(grepl("-|7SK", g)){
      g <- paste0("`",g,"`")
    }
    # snp-gene
    snp.cis_gene <- lm(paste0(g, "~",snp), data=d)
    snp.cis_gene <- coefficients(snp.cis_gene)[snp]

    # snp-trans_assoc and cis_gene-trans_assoc
    for(ta in trans_assoc) {
      if(grepl("-|7SK", ta)){
        ta <- paste0("`",ta,"`")
      }
      snp.trans_assoc <- lm(paste0(ta,"~",snp), data=d)
      snp.trans_assoc <- coefficients(snp.trans_assoc)[snp]

      cis_gene.trans_assoc <- lm(paste0(ta, "~", g), data=d)
      cis_gene.trans_assoc <- coefficients(cis_gene.trans_assoc)[g]

      # snp-cpg hat
      snp.trans_assoc.hat <- snp.cis_gene * cis_gene.trans_assoc

      # save result
      betas <- rbind(betas,
                     c(snp, g, ta, snp.trans_assoc, snp.cis_gene,
                       cis_gene.trans_assoc, snp.trans_assoc.hat))
    }
  }
  colnames(betas) <- c("snp", "cis_gene", "trans_assoc", "snp.trans_assoc",
                       "snp.cis_gene", "cis_gene.trans_assoc", "snp.trans_assoc.hat")
  betas <- as.data.frame(betas, stringsAsFactors=F)
  betas$snp.trans_assoc <- as.numeric(betas$snp.trans_assoc)
  betas$snp.cis_gene <- as.numeric(betas$snp.cis_gene)
  betas$cis_gene.trans_assoc <- as.numeric(betas$cis_gene.trans_assoc)
  betas$snp.trans_assoc.hat <- as.numeric(betas$snp.trans_assoc.hat)

  # plot the betas for each gene
  gp <- ggplot(data=betas, aes(x=snp.trans_assoc, y=snp.trans_assoc.hat)) +
    geom_hline(yintercept=0, colour="grey") +
    geom_vline(xintercept=0, colour="grey") +
    facet_wrap(~ cis_gene, ncol=3) +
    geom_point() +
    geom_smooth(method="lm") +
    ylab(expression(hat(beta)[c])) +
    xlab(expression(beta[c])) +
    geom_abline(slope=1, colour="orange")
  ggsave(fout_plot, plot=gp, width=12, height=12)

  # also save the list of betas
  write.table(file=fout_table,
              betas, col.names=T, sep="\t", quote=F,
              row.names=F)

  # get the correlation results for the estimated against the
  # observed betas
  result <- lapply(cis_genes, function(g){
    if(grepl("-|7SK", g)){
      g <- paste0("`",g,"`")
    }

    d <- betas[betas$cis_gene == g,]
    d1 <- d[,"snp.trans_assoc"]
    d2 <- d[,"snp.trans_assoc.hat"]
    fit <- lm(d2~d1)
    r <- cor.test(d1, d2)
    cor <- r$estimate
    pv <- r$p.value

    # also construct a contingency table where we only look at the signs
    # of the individual beta coefficients for each CpG and test the association
    tab <- matrix(0,ncol=2,nrow=2)
    for(i in 1:nrow(d)) {
      obs <- d[i,"snp.trans_assoc"]
      pred <- d[i,"snp.trans_assoc.hat"]
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
  names(result) <- cis_genes
  return(result)
}

#' -----------------------------------------------------------------------------
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
#' -----------------------------------------------------------------------------
mediation.summary <- function(med, s, s.selected, cutoff) {
  # mediation for only ggm selected snp genes
  med.selected <- med[s.selected]
  med.selected.pvals <- unlist(lapply(med.selected, "[[", "pval_fisher"))
  med.selected.sign <- med.selected[med.selected.pvals<cutoff]

  # mediation for only not selected snp genes
  med.notselected <- med[setdiff(s,s.selected)]
  med.notselected.pvals <- unlist(lapply(med.notselected, "[[", "pval_fisher"))
  med.notselected.sign <- med.notselected[med.notselected.pvals<cutoff]

  # get the best mediating gene based on the correlation
  best <- med[which.max(abs(unlist(lapply(med, "[[", "odds_ratio"))))]
  best_mediator <- names(best)
  best_odds_ratio <- unname(best[[1]]["odds_ratio"])

  # TODO if we choose to do a test for comparing
  # selected/not selected, check assumptions first

  # min mediation pv for all not selected snpgenes
  med.pv <- min(med.notselected.pvals)
  # max mediation pv for selected snpgenes
  med.pv.selected <- max(med.selected.pvals)

  # number of significantly mediating genes
  med_total <- length(which(c(med.notselected.pvals,med.selected.pvals)<cutoff))

  cat("Mediation result summary (min/max of pvals):\n")
  cat("not selected\tselected\n")
  cat(med.pv, "\t", med.pv.selected, "\n")

  # compare the two pvalues
  d <- log10(med.pv/med.pv.selected)
  cat("Difference (log10(ns/s)): ", d, "\n")

  # return 'nice' result
  return(c(mediation_total=med_total,
           mediation_best_gene=best_mediator,
           mediation_best_odds_ratio=best_odds_ratio,
           mediation_notselected_significant=length(med.notselected.sign),
           mediation_selected_significant=length(med.selected.sign),
           mediation_notselected_significant.list=paste(names(med.notselected.sign), collapse=";"),
           mediation_selected_significant.list=paste(names(med.selected.sign), collapse=";"),
           mediation_notselected.list=paste(names(med.notselected), collapse=";"),
           mediation_selected.list=paste(names(med.selected), collapse=";"),
           mediation_notselected_pvals=paste(format(med.notselected.pvals, digits=3), collapse=";"),
           mediation_selected_pvals=paste(format(med.selected.pvals, digits=3), collapse=";"),
           mediation_min_pval_notselected=med.pv,
           mediation_max_pval_selected=med.pv.selected,
           log10_mediation_NSoverS_ratio=d))
}

#' -----------------------------------------------------------------------------
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
#' -----------------------------------------------------------------------------
compare_mediation_results <- function(sentinel,
                                      mediation,
                                      mediation2,
                                      sgenes_selected,
                                      cutoff,
                                      fout=NULL) {
  require(ggplot2)
  require(reshape2)
  library(cowplot)

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

  df <- df[complete.cases(df),,drop=F]
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
  df <- reshape2::melt(df, measure.vars=1:2, variable.name="dataset")

  # add plotting x-lab proxy
  df = transform(df,
                 dvariable = as.numeric(dataset) + runif(length(dataset),
                                                          -0.2,
                                                          0.2))

  # create summary plot
  #gp <- ggplot(data=df, aes(y=value, x=dataset))
  #gp <- gp + geom_violin(draw_quantiles=.5) +
  #  geom_point(aes(color=significant, x=dvariable, shape=model_selected), size=3) +
  #  scale_color_brewer(palette = "Set2") +
  #  geom_line(aes(group=gene, x=dvariable), linetype="dotted")+
  #  geom_text(aes(label=gene), size=2.5, nudge_x=-0.4,
  #            data=subset(df, dataset=="correlation_validation")) +
  #  ggtitle(paste0("Mediation correlation values for ", sentinel, "."),
  #          "Shown are the distributions of correlation values over all SNP genes,
#validation dataset VS the original dataset (on which models were calculated)") +
  #  theme(plot.title = element_text(hjust = 0)) + background_grid(major = "xy")

  # save to file
  #ggsave(file=fout,
  #       plot=gp)

  # return the correlation of mediation correlations between the two datasets
  # and the fraction of significant mediations in validation set also
  # significant in original dataset
  corr <- cor(med_correlations, med_correlations2)
  df <- subset(df, dataset=="correlation_validation")
  # fraction of all SNP genes significant in both datasets
  frac <- length(which(df$significant=="both"))/nrow(df)
  # fraction of validation significant SNP genes significant in original dataset
  frac2 <- length(which(df$significant=="both" | df$significant == "original")) /
    length(which(df$significant=="both" | df$significant=="validation"))

  return(c(mediation_cross_cohort_correlation=corr,
           mediation_cross_cohort_fraction=frac,
           mediation_cross_cohort_fraction_validation_significant=frac2))
}
