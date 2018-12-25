# ------------------------------------------------------------------------------
#' Collects methylation, expression and genotype data for a either KORA or
#' LOLIPOP cohort based on a ranges object
#'
#' @autho Johann Hawe
#'
#' @date 02/06/2017
#'
#' @export
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


# ------------------------------------------------------------------------------
print("Load libraries, source scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# Define methods
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#' Adjust the data collected for a specific sentinel SNP
#'
#' Adjusts all data for a given sentinel SNP, based on the previously
#' identified ranges, genes etc done in \code{collect.ranges}.
#' Needs methylation and expression data available
#' ('meth' and 'expr' variables, respectively).
#'
#' @param sentinel the id of the sentinel SNP to be processed
#' @param ranges The collected ranges for the sentinel SNP
#' @param data A Matrix containing the data retrieved from the ranges object
#'
#' @value NULL
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
adjust_data <- function(sentinel, ranges, data, geno, fccosmo, fceqtl) {

  # collect gene symbols to get summarized probe levels
  symbols <- unique(c(ranges$cpg_genes$SYMBOL,
                      ranges$snp_genes$SYMBOL,
                      ranges$spath$SYMBOL,
                      ranges$tfs$SYMBOL,
                      ranges$trans_genes$SYMBOL))
  symbols <- symbols[!sapply(symbols, is.na)]

  # retrieve the genotype data
  s <- data[,sentinel,drop=F]

  # correct for covariate variation
  g_resid <- rm_covariate_effects(data[,!grepl("^rs|^cg",
                                               colnames(data))], "expr")

  # ----------------------------------------------------------------------------
  # meQTL seed sepcific processing
  # ----------------------------------------------------------------------------
  if(ranges$seed == "meqtl") {
    # remember cpg genes' probe ids -> used for adjusting cis eqtl effects
    cpg_gene_probes <- unique(unlist(ranges$cpg_genes$ids))

    # process gene data
    print("Loading eQTLs.")
    eqtls <- load_eqtls(fceqtl, cpg_gene_probes)

    # get genotype data for eqtls and meqtls
    ids <- unique(eqtls$snp_id)
    if(!is.null(ids)) {
      gen <- geno[,colnames(geno) %in% ids, drop=F]
      gen <- gen[data$geno_ids,,drop=F]

      # get rid of cis-eqtl effects
      print("Adjusting for cis eqtls.")
      g_resid <- adjust_cis_eqtls(g_resid, eqtls, gen)
    }

    # process the CpG data
    c_resid <- rm_covariate_effects(data[,!grepl("^rs",
                                                 colnames(data))], "meth")

    print("Loading meQTLs.")
    meqtls <- load_meqtls(fccosmo, colnames(c_resid))
    ids <- meqtls$snp_id
    if(!is.null(ids)) {
      gen <- geno[, colnames(geno) %in% ids, drop=F]
      gen <- gen[data$geno_ids,,drop=F]
      # get rid of cis-eqtl effects
      print("Adjusting for cis meqtls.")
      c_resid <- adjust_cis_meqtls(c_resid, meqtls, gen)
    }
  }

  print("Summarizing probe levels.")
  # summarize to gene level estimates from expression probes
  g_resid <- summarize(g_resid, symbols = symbols)

  # create complete matrix, containing all the information
  if(ranges$seed == "meqtl") {
    data <- cbind.data.frame(g_resid, c_resid, s,
                             stringsAsFactors=F)
  } else {
    data <- cbind.data.frame(g_resid, s,
                             stringsAsFactors=F)
  }
  return(data)
}

# ------------------------------------------------------------------------------
#' Gets residuals based on linear models from a data matrix
#'
#' Calculates the residual matrix for a given matrix considering available
#' covariates.
#' Uses linear model for methylation data and linear mixed
#' model with plate as random effect for the expression data.
#'
#' @param data the matrix for which to calculate the covariates
#' @param data.type the type of data in the matrix: either "meth" or "expr".
#' Depending on the type different formulas are used for the linear model.
#' @param cols The col.names over which to iterate in the dataframe to calculate
#' the residuals on (e.g. probe.ides, gene.names,..)
#'
#' @return A matrix  where in the colums are the measured entities (e.g.probes)
#' and in the rows are the samples, containing the calculated residuals
#'
# ------------------------------------------------------------------------------
rm_covariate_effects <- function(data, data.type, cols=NULL) {

  if(data.type=="meth") {
    if(is.null(cols)) {
      cols <- colnames(data)[grepl("^cg", colnames(data))]
    }
    # define design and response matrix
    design <- data[,c("CD4T", "CD8T",
		      "NK", "Bcell", "Mono",
		      paste("PC", paste(1:20, "cp", sep="_"), sep=""))]
    resp <- data[,cols]
  } else if(data.type == "expr") {
    if(is.null(cols)){
      cols <- colnames(data)[grepl("^ILMN_", colnames(data))]
    }
    # define design and response matrix
    design <- data[,c("age", "sex", "RIN", "batch1", "batch2")]
    resp <- data[,cols]
  } else {
    stop("Data type not supported for residual calculation.")
  }

  resp <- data.matrix(resp)
  design <- data.matrix(design)

  # calculate lm and get residuals
  residuals <- resid(lm(resp~design, na.action=na.exclude))

  return(residuals)
}

# ------------------------------------------------------------------------------
#' Takes a gene expression matrix and adjusts the genes' expression
#' for cis-eQTLs
#'
#' @param expr Gene expression matrix with all the genes in the column.
#' Column names need to be illumina probe ids.
#' @param eqtls List of eqtls for which to adjust the data
#' @param geno_data Available genotype data
#'
#' @return The gene expression matrix adjusted for cis-eQTLs
#'
#' @author  Johann Hawe
#'
# ------------------------------------------------------------------------------
adjust_cis_eqtls <- function(expr, eqtls, geno_data){

  # for all cpg-genes, adjust for potential cis-eqtls (i.e. snp-effects)
  toUse <- which(eqtls$probe_id %in% colnames(expr))
  if(length(toUse) == 0){
    return(expr)
  }

  eqtls <- eqtls[toUse, , drop=F]
  # for each gene, check whether we have an associated cis-eqtl
  temp <- lapply(colnames(expr), function(x) {
    # indices of ciseqtls
    idxs <- which((eqtls$probe_id == x) & eqtls$snp_id %in% colnames(geno_data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    snpids <- eqtls[idxs,,drop=F]$snp_id
    X <- cbind.data.frame(expr[,x,drop=F], geno_data[,snpids,drop=F])
    model <- lm(as.formula(paste0("`", x, "` ~",
                                  paste0(snpids,collapse="+"))),
                data =  X,
                na.action=na.exclude)

    r <- resid(model)
    return(r)
  })

  names(temp) <- colnames(expr)

  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      expr[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }

  print(paste0("Adjusted expression value of ", cnt, " genes."))
  return(expr)
}

# ------------------------------------------------------------------------------
#' Takes a methylation matrix and adjusts the cpgs' expression for cis-meQTLs
#' using currently the BONDER meqtl results.
#'
#' @param betas Methylation beta matrix with all the cpgs in the columns.
#' @param meqtls The identified cis meQTLs
#' @param geno_data The available genotype data
#'
#' @return The methylation matrix adjusted for cis-meQTLs
#'
#' @author  Johann Hawe
#'
# ------------------------------------------------------------------------------
adjust_cis_meqtls <- function(betas, meqtls, geno_data){
  # for all cpg-genes, adjust for potential cis-eqtls (i.e. snp-effects)
  toUse <- which(meqtls$cpg_id %in% colnames(betas))
  if(length(toUse) == 0){
    return(betas)
  }

  meqtls <- meqtls[toUse,]

  # for each gene, check whether we have an associated cis-meqtl
  temp <- lapply(colnames(betas), function(x) {
    # indices of cis-meqtls
    idxs <- which((meqtls$cpg_id == x) & meqtls$snp_id %in% colnames(geno_data))
    if(length(idxs) < 1) {
      return(NULL)
    }
    snpids <- meqtls[idxs,,drop=F]$snp_id
    X <- cbind.data.frame(betas[,x,drop=F],
                          geno_data[,snpids,drop=F])

    # drop data with no variance
    X <- X[,unlist(apply(X,2,function(x){ var(x,na.rm = T)!=0 })), drop=F]
    if(ncol(X)<1 | !x %in% colnames(X)){
      return(NULL)
    }
    model <- lm(as.formula(paste0("`", x, "` ~",
                                  paste0(snpids[snpids %in% colnames(X)],
                                         collapse="+"))),
                data =  X,
                na.action=na.exclude)

    r <- resid(model)
    return(r)
  })

  names(temp) <- colnames(betas)

  # set results to our data matrix
  cnt <- 0
  for(i in 1:length(temp)){
    x <- names(temp)[i]
    if(!is.null(temp[[i]])) {
      betas[,x] <- temp[[i]]
      cnt <- cnt + 1
    }
  }

  print(paste0("Adjusted methylation beta of ", cnt, " cpgs."))
  return(betas)
}

# ------------------------------------------------------------------------------
#' Load cis-eqtls identified in KORA
#'
#' Loads a list of cis-eqtls identified previously in the KORA F4 cohort.
#' (Supplement Table S2 of Schramm et al.).
#'
#' @param fceqtl Path to file containing cis eqtl
#' @param expr.probes Optional. Restrict the eqtls to be loaded to the
#' this list of expression probes obly
#'
#' @author Johann Hawe
#'
#' @return GRanges objects containing the eqtl information
#'
# ------------------------------------------------------------------------------
load_eqtls <- function(fceqtl, expr.probes=NULL) {
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

  # check whether we got any hits
  if(length(genes) < 1) {
    return(NULL)
  }

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

# ------------------------------------------------------------------------------
#' Load cis meqtls identified in LOLIPOP/KORA study
#'
#' @param fccosmo Path to the file containing all available cis-meQTL
#' @param cpg_probes Optional. Restrict the list of meqtls on those
#' cpg-probes only.
#'
#' @author Johann Hawe
#'
#' @return GRanges object reflecting all found cis-meQTLs
#'
# ------------------------------------------------------------------------------
load_meqtls <- function(fccosmo, cpg_probes=NULL) {

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

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# input
franges <- snakemake@input[["ranges"]]
fkora_data <- snakemake@input[["kora"]]
flolipop_data <- snakemake@input[["lolipop"]]
fccosmo <- snakemake@input[["ccosmo"]]
fceqtl <- snakemake@input[["ceqtl"]]

# params
cohort <- snakemake@wildcards$cohort
sentinel <- snakemake@wildcards$sentinel
seed <- snakemake@wildcards$seed

# output
fout <- snakemake@output[[1]]
fout_raw <- snakemake@output[[2]]

print(paste0("Sentinel is ", sentinel, "."))

# ------------------------------------------------------------------------------
print("Prepare probe and snp ids.")
# ------------------------------------------------------------------------------
ranges <- readRDS(franges)
expr_probes = c(unlist(ranges$snp_genes$ids),
                unlist(ranges$cpg_genes$ids),
                unlist(ranges$tfs$ids),
                unlist(ranges$spath$ids))
expr_probes <- unique(expr_probes)
meth_probes <- names(ranges$cpgs)

# ------------------------------------------------------------------------------
print("Loading cohort data.")
# ------------------------------------------------------------------------------
if("kora" %in% cohort) {
  load(fkora_data)
} else if ("lolipop" %in% cohort) {
  load(flolipop_data)
} else {
  stop("Cohort not supported.")
}

# ------------------------------------------------------------------------------
print("Create merged data frame.")
# ------------------------------------------------------------------------------
data <- cbind.data.frame(covars,
                         expr[,colnames(expr) %in% expr_probes, drop=F],
                         meth[,colnames(meth) %in% meth_probes, drop=F],
                         geno[,colnames(geno) %in% sentinel, drop=F],
                         geno_ids=rownames(geno),
                         stringsAsFactors=F)

# ------------------------------------------------------------------------------
print("Save raw data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout_raw, data)

# check which probes we didnt have available
print("Unavailable probes:")
print(expr_probes[!expr_probes %in% colnames(expr)])
print(meth_probes[!meth_probes %in% colnames(meth)])

# remove not needed remaining data frames
rm(expr, meth, covars)
gc()

print(paste0("Data dimensions: ", paste(dim(data), collapse=",")))

# ------------------------------------------------------------------------------
print("Removing covariate effects from raw data.")
# ------------------------------------------------------------------------------
data <- adjust_data(sentinel, ranges, data, geno, fccosmo, fceqtl)

# ------------------------------------------------------------------------------
print("Saving adjusted data.")
# ------------------------------------------------------------------------------
saveRDS(file=fout, data)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
