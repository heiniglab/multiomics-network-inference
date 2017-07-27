#' Performs validation of given ggm results.
#' 
#' @author Johann Hawe
#' 
#' @version 20170720
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
#' @version 20170719
#' 
mediation <- function(data, snp, genes, cpgs, plot=T) {
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
    r <- lm(paste0(g,"~",snp), data=d)
    coefficients(r)[snp]
  })
  names(snp.gene) <- genes
  
  # snp-cpg
  snp.cpg <- lapply(cpgs, function(c){
    r <- lm(paste0(c,"~",s), data=d)
    coefficients(r)[snp]
  })
  names(snp.cpg) <- cpgs
  
  betas <- lapply(genes, function(g){
    # gene-cpg
    temp <- lapply(cpgs, function(c){
      r <- lm(paste0(c,"~",g), data=d)
      # merge with the other betas (for each cpg and gene combi, we have 3 betas)
      r <- c(snp.cpg[[c]], snp.gene[[g]], coefficients(r)[g])
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

# prepare data
load("results/rs9859077.test.fit.RData")
data <- data.test
snp <- names(ranges$sentinel.range)
data[,snp] <- as.integer(as.character(data[,snp]))
cpgs <- names(ranges$cpgs)
genes <- ranges$snp.genes$SYMBOL
genes <- genes[genes %in% colnames(data)]

gn <- nodes(graph)
selected.genes <- gn[unlist(nodeData(graph,gn,"snp.gene"))]

par(mfrow=c(4,round(length(genes)/4+0.5)), mar=rep(4.5,4))
pdf("senp7.mediation.pdf")
m <- mediation(data, snp, genes, cpgs)
dev.off()
