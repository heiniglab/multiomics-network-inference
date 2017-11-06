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
      plot(x = d1,y = d2,
           main=paste0(g),
           ylab=expression(hat(beta)[c]),
           xlab=expression(beta[c]),
           pch="+")
      abline(fit, col="red")
      text(x=-0.01, y=max(d1), labels = paste0("cor=", format(cor, digits=2), 
                             "\np.value=", format(pv, digits=2)))
    }
    c(correlation=unname(cor),pvalue=pv)
  })
  names(result) <- genes
  result
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

#' Validates the given ggm fit based on the found Gene~Gene links
#' 
#' @param ggm.fit The ggm.fit which is to be validated
#' @param ranges The original list of granges related to the entities in the 
#' fitted data matrix
#' 
#' @return 
#' 
#' @author Johann Hawe
#' 
validate.genes <- function(ggm.fit, ranges){
  
}

sigma.trace <- function(ggm.fit) {
  
}
