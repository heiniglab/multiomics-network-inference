#' Collects methylation, expression and genotype data for a either KORA or LOLIPOP
#' cohort. 
#' 
#' @autho Johann Hawe
#' 
#' @date 02/06/2017
#' 
#' @export
#' 

library(GenomicRanges)
source("scripts/lib.R")

#' Adjust the data collected for a specific sentinel SNP
#'
#' Adjusts all data for a given sentinel SNP, based on the previously
#' identified ranges, genes etc done in \code{collect.ranges}. 
#' Needs methylation and expression data available
#' ('meth' and 'expr' variables, respectively). 
#'
#' @param sentinel the id of the sentinel SNP to be processed (needs to be in the list "sentinels")
#' @param ranges The collected ranges for the sentinel SNP
#' @param data A Matrix containing the data retrieved from the ranges object
#'
#' @value NULL
#' 
#' @author Johann Hawe
#'
adjust.data <- function(sentinel, ranges, data, geno, fccosmo, fceqtl) {
  
  # collect gene symbols to get summarized probe levels
  symbols <- unique(c(ranges$cpg_genes$SYMBOL, 
                      ranges$snp_genes$SYMBOL,
                      ranges$spath$SYMBOL,
                      ranges$tfs$SYMBOL))
  symbols <- symbols[!sapply(symbols, is.na)]
  
  # retrieve the genotype data
  s <- data[,sentinel,drop=F]
  
  # correct for covariate variation
  g_resid <- get.residuals(data[,!grepl("^rs|^cg", colnames(data))],"expr")
  c_resid <- get.residuals(data[,!grepl("^rs", colnames(data))],"meth")
  
  cat("Loading QTLs.\n")
  eqtls <- load.eqtls(fceqtl, colnames(g_resid))
  meqtls <- load.meqtls(fccosmo, colnames(c_resid))
  
  # get genotype data for eqtls and meqtls
  ids <- unique(eqtls$snp_id)
  if(!is.null(ids)) {
    gen <- geno[,colnames(geno) %in% ids, drop=F]
    gen<- gen[data$geno_ids,,drop=F]
    
    # get rid of cis-eqtl effects
    cat("Adjusting for cis eqtls.\n")
    g_resid <- adjust.cis.eqtls(g_resid, eqtls, gen)
  }
  ids <- meqtls$snp_id
  if(!is.null(ids)) {
    gen <- geno[, colnames(geno) %in% ids, drop=F]
    gen<- gen[data$geno_ids,,drop=F]
    # get rid of cis-eqtl effects
    cat("Adjusting for cis meqtls.\n")
    c_resid <- adjust.cis.meqtls(c_resid, meqtls, geno)
  }
  
  cat("Summarizing probe levels.\n")
  # summarize to gene level estimates from expression probes
  g_resid <- summarize(g_resid, symbols = symbols)
  
  # create complete matrix, containing all the information
  data <- cbind.data.frame(g_resid,
                           c_resid,
                           s,
                           stringsAsFactors=F)
  return(data)
}

#' 
#' Method to get gene level estimates from expression probes
#' (=rowMeans on respective probe.ids/gene)
#' 
#' @param m The expression matrix with probe ids in columns
#' @param symbols The gene symbols to be used for the summarization. Only probe ids
#' belonging to a gene symbol in this list a used. Rest of probes is dropped.
#' 
#' @value Matrix of summarized probe values (per gene probe levels)
#' 
#' @author Johann Hawe
#' 
#' @date 02/06/2017
#' 
#' @export
#' 
summarize <- function(m, symbols){
  # prepare annotation
  library(illuminaHumanv3.db)
  annot <- as.list(illuminaHumanv3ALIAS2PROBE)
  annot <- annot[!is.na(annot)]
  n <- lapply(symbols, function(sym) {
    gid <- probes.from.symbols(sym, annot=annot)
    # for LST1 symbol (and possibly others) there is a probe which should likely
    # be only appointed to the NFKBIL1 gene (according to UCSC genome browser). 
    # Might be an error in annotation. Remove probe for LST1 manually for now
    if(sym=="LST1"){
      gid <- gid[which(!(gid %in% "ILMN_2046344"))]
    }
    if(is.null(gid)){
      return(NULL)
    }
    if(length(gid)>0){
      if(!all(gid %in% colnames(m))){
        gid <- gid[which(gid %in% colnames(m))]
        if(length(gid)<1){
          warning("None of the probeids for ", 
                  sym, 
                  " where available in expression data.")
          return(NULL)
        } else {
          warning("Some probes where not available in expression data.")
        }
      }
      
      g <- m[,gid,drop=F]
      s <- as.matrix(rowMeans(g), ncol=1)
      colnames(s) <- sym
      return(s)
    }
    return(NULL)
  })
  names(n) <- symbols
  n <- n[!sapply(n,is.null)]
  result <- matrix(unlist(n), ncol=length(n), byrow=F)
  colnames(result) <- names(n)
  rownames(result) <- rownames(m)
  return(result)
}

#' Gets residuals based on linear models from a data matrix
#'
#' Calculates the residual matrix for a given matrix considering available covariates. 
#' Uses linear model for methylation data and linear mixed
#' model with plate as random effect for the expression data.
#'
#' @param data the matrix for which to calculate the covariates
#' @param data.type the type of data in the matrix: either "meth" or "expr". Depending on the type
#' different formulas are used for the linear model.
#' @param col.names The col.names over which to iterate in the dataframe to calculate
#' the residuals on (e.g. probe.ides, gene.names,..)
#' 
#' @return A matrix  where in the colums are the measured entities (e.g.probes) 
#' and in the rows are the samples, containing the calculated residuals
#'
get.residuals <- function(data, data.type, col.names=NULL) {
  
  cols <- col.names
  
  if(data.type=="meth") {
    if(is.null(cols)) {
      cols <- colnames(data)[grepl("^cg", colnames(data))]
    }
    res <- lapply(cols, function(n) {
      fm <- as.formula(paste0(n,"~",
                              paste0("1+CD4T+CD8T+NK+Bcell+Mono+",
                                     paste0("PC", 
                                            paste(1:20, "cp", sep="_"), 
                                            collapse="+"))))
      return(lm(fm, data=data, na.action=na.exclude))
    })
  } else if(data.type == "expr") {
    if(is.null(cols)){
      cols <- colnames(data)[grepl("^ILMN_", colnames(data))]
    }
    res <- lapply(cols, function(n) {
      fm <- as.formula(paste0(n, "~",
                              "1+age+sex+RIN+batch1+batch2"))
      return(lm(fm,data=data, na.action=na.exclude))
    })
  } else {
    stop("Data type not supported for residual calculation.")
  }
  # build the full residual matrix from model results
  residual.mat <- matrix(data=unlist(lapply(res, resid)), nrow=nrow(data))
  colnames(residual.mat) <- cols
  rownames(residual.mat) <- rownames(data)
  
  return(residual.mat)
}

#' Takes a gene expression matrix and adjusts the genes' expression for cis-eQTLs
#' using gtex eqtl results.
#'
#' @param g Gene expression matrix with all the genes in the column. Column names
#' need to be gene symbols.
#' @param eqtls List of eqtls for which to adjust the data
#'
#' @return The gene expression matrix adjusted for cis-eQTLs
#'
#' @author  Johann Hawe
#'
adjust.cis.eqtls <- function(g, eqtls, geno_data){
  
  # for all cpg-genes, adjust for potential cis-eqtls (i.e. snp-effects)
  toUse <- which(eqtls$gene %in% colnames(g))
  if(length(toUse) == 0){
    return(g)
  }
  
  eqtls <- eqtls[toUse, , drop=F]
  # for each gene, check whether we have an associated cis-eqtl
  temp <- lapply(colnames(g), function(x) {
    # indices of ciseqtls
    idxs <- which((eqtls$gene == x) & eqtls$snp_id %in% colnames(geno_data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    snpids <- eqtls[idxs,,drop=F]$snp_id
    X <- cbind.data.frame(g[,x,drop=F], geno_data[,snpids,drop=F])
    model <- lm(as.formula(paste0("`", x, "` ~",
                                  paste0(snpids,collapse="+"))), 
                data =  X,
                na.action=na.exclude)
    
    r <- resid(model)
    return(r)
  })
  
  names(temp) <- colnames(g)
  
  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      g[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }
  
  cat("Adjusted expression value of", cnt, "genes.\n")
  return(g)
}

#' Takes a methylation matrix and adjusts the cpgs' expression for cis-meQTLs
#' using currently the BONDER meqtl results.
#'
#' @param c Methylation beta matrix with all the cpgs in the columns.
#'
#' @return The methylation matrix adjusted for cis-meQTLs
#'
#' @author  Johann Hawe
#'
adjust.cis.meqtls <- function(c, meqtls, geno_data){
  # for all cpg-genes, adjust for potential cis-eqtls (i.e. snp-effects)
  toUse <- which(meqtls$cpg_id %in% colnames(c))
  if(length(toUse) == 0){
    return(c)
  }
  
  meqtls <- meqtls[toUse,]
  
  # for each gene, check whether we have an associated cis-meqtl
  temp <- lapply(colnames(c), function(x) {
    # indices of cis-meqtls
    idxs <- which((meqtls$cpg_id == x) & meqtls$snp_id %in% colnames(geno_data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    snpids <- meqtls[idxs,,drop=F]$snp_id
    X <- cbind.data.frame(c[,x,drop=F], geno_data[,snpids,drop=F])
    # drop data where all values are the same
    X <- X[,unlist(apply(X,2,function(x){ var(x,na.rm = T)!=0 })), drop=F]
    if(ncol(X)<1 | !x %in% colnames(X)){
      return(NULL)
    }
    model <- lm(as.formula(paste0("`", x, "` ~",
                                  paste0(snpids[snpids %in% colnames(X)],collapse="+"))), 
                data =  X,
                na.action=na.exclude)
    
    r <- resid(model)
    return(r)
  })
  
  names(temp) <- colnames(c)
  
  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      c[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }
  
  cat("Adjusted methylation beta of", cnt, "cpgs\n")
  return(c)
}

#' Load cis-eqtls identified in KORA
#'
#' Loads a list of cis-eqtls identified previously in the KORA F4 cohort.
#' (Supplement Table S2 of Schramm et al.).
#'
#' @param expr.probes Optional. Restrict the eqtls to be loaded to the this list of expression
#' probes obly
#' 
#' @author Johann Hawe
#' 
#' @return GRanges objects containing the eqtl information
#'
load.eqtls <- function(fceqtl, expr.probes=NULL) {
  # columns: top SNP KORA F4;Chr SNP;minor allele KORA F4;MAF KORA F4;Probe_Id;Gene;Chr Gene; etc...
  eqtls <- read.csv(fceqtl,
                    sep=";", header=T, stringsAsFactors=F)
  
  sids <- eqtls$top.SNP.KORA.F4
  keep <- which(rep(T,length(sids)))
  if(!is.null(expr.probes)){
    keep <- which(eqtls$Probe_Id %in% expr.probes)
  }
  
  sids <- sids[keep]
  genes <- eqtls$Gene[keep]
  pids <- eqtls$Probe_Id[keep]
  
  # try to map some additional symbols
  temp <- lapply(1:length(genes), function(i) {
    if(genes[i] == "") {
      sym <- symbols.from.probeids(pids[i])
      genes[i] <<- sym
    }
  })
  
  # drop duplicates (since the analysis was probebased and we are only on the
  # per gene level)
  keep <- which(!duplicated(paste0(genes, sids)))
  eqtls <- cbind.data.frame(snp_id=sids[keep],
                            stringsAsFactors=F)
  rownames(eqtls) <- paste0(genes[keep], sids[keep])
  eqtls$probe_id <- pids[keep]
  eqtls$gene <- genes[keep]
  return(eqtls)
}

#' Load cis meqtls identified in LOLIPOP/KORA study
#' 
#' @param cpg_probes Optional. Restrict the list of meqtls on those cpg-probes only.
#' 
#' @author Johann Hawe
#' 
#' @return GRanges object reflecting all found cis-meQTLs
#' 
load.meqtls <- function(fccosmo, cpg_probes=NULL) {
  
  cosmo <- readRDS(fccosmo)
  
  if(is.null(cpg_probes)){
    cosmosub <- cosmo
  } else {
    cosmosub <- cosmo[(cosmo$cpg %in% cpg_probes),,drop=F]
  }
  rm(cosmo)
  gc()
  
  if(nrow(cosmosub) > 0){
    # load snp ranges (columns: chr pos id)
    meqtl <- with(cosmosub, 
                  cbind.data.frame(snp_id=as.character(snp), 
                                   cpg_id=as.character(cpg),
                                   stringsAsFactors=F))
    rownames(meqtl) <- paste0(meqtl$snp_id, meqtl$cpg_id)
    
    return(meqtl)
  } else {
    return(NULL)
  }
}

# start processing -------------------------------------------------------------

# get params -------------------------------------------------------------------
franges <- snakemake@input[["ranges"]]
fkora_data <- snakemake@input[["kora"]]
flolipop_data <- snakemake@input[["lolipop"]]
fccosmo <- snakemake@input[["ccosmo"]]
fceqtl <- snakemake@input[["ceqtl"]]
cohort <- snakemake@params$cohort
sentinel <- snakemake@params$sentinel
fout <- snakemake@output[[1]]

print(paste0("Sentinel is ", sentinel, "."))

# prepare probe and snp ids ----------------------------------------------------
ranges <- readRDS(franges)
expr_probes = c(unlist(ranges$snp_genes$ids), 
                unlist(ranges$cpg_genes$ids))
if("tfs" %in% names(ranges)){
  expr_probes <- c(expr_probes, unlist(ranges$tfs$ids))
}
if("spath" %in% names(ranges)){
  expr_probes <- c(expr_probes, unlist(ranges$spath$ids))
}
expr_probes <- unique(expr_probes)
meth_probes <- names(ranges$cpgs)
snp_ids <- sentinel

print("Loading cohort data.")
if("kora" %in% cohort) {
  load(fkora_data)   
} else if ("lolipop" %in% cohort) {
  load(flolipop_data)
} else {
  stop("Cohort not supported.")
}

# create merged data frame -----------------------------------------------------

data <- cbind.data.frame(covars, 
                         expr[,colnames(expr) %in% expr_probes, drop=F], 
                         meth[,colnames(meth) %in% meth_probes, drop=F], 
                         geno[,colnames(geno) %in% snp_ids, drop=F], 
                         geno_ids=rownames(geno),
                         stringsAsFactors=F)
print(colnames(data))
# remove not needed remaining data frames
rm(expr, meth, covars)
gc()

print(paste0("Data dimensions: ", paste(dim(data), collapse=",")))
data <- adjust.data(sentinel, ranges, data, geno, fccosmo, fceqtl)

saveRDS(file=fout, data)
