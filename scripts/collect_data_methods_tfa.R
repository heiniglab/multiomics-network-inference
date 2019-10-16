#' -----------------------------------------------------------------------------
#' Convenience methods for the collect data script
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Jun  4 17:12:58 2019
#' -----------------------------------------------------------------------------

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
adjust_data <- function(sentinel, ranges, genes, data, geno, fccosmo, fceqtl) {

  # retrieve the genotype data
  s <- data[,sentinel,drop=F]
  g <- data[,genes,drop=F]

  # ----------------------------------------------------------------------------
  # meQTL seed specific processing
  # ----------------------------------------------------------------------------
  if(ranges$seed == "meqtl") {
    # remember cpg genes' probe ids -> used for adjusting cis eqtl effects
    cpg_genes <- unique(unlist(ranges$cpg_genes$SYMBOL))
    cpg_genes <- cpg_genes[cpg_genes %in% genes]

    # process gene data
    print("Loading eQTLs.")
    eqtls <- load_eqtls(fceqtl, cpg_genes)

    # get genotype data for eqtls and meqtls
    ids <- unique(eqtls$snp_id)
    if(!is.null(ids)) {
      gen <- geno[,colnames(geno) %in% ids, drop=F]
      gen <- gen[data$geno_ids,,drop=F]

      # get rid of cis-eqtl effects
      print("Adjusting for cis eqtls.")
      g <- adjust_cis_eqtls(g, eqtls, gen, cpg_genes)
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

  # create complete matrix, containing all the information
  if(ranges$seed == "meqtl") {
    # scale and center
    resid <- scale(cbind.data.frame(g, c_resid))
    data <- cbind.data.frame(resid, s,
                             stringsAsFactors=F)
  } else {
    resid <- scale(g)
    data <- cbind.data.frame(resid, s,
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

  residuals <- lapply(colnames(resp), function(v) {
    resid(lm(resp[,v]~design, na.action=na.exclude))
  })

  # calculate lm and get residuals
  # NOTE: we tried resp ~ design as a single call before, however, due to
  # potential NAs in the data, this lead to ALL non complete.cases being removed
  # before doing the fits, which we dont want here
  residuals <- do.call(cbind, residuals)
  colnames(residuals) <- colnames(resp)

  return(residuals)
}

# ------------------------------------------------------------------------------
#' Takes a gene expression matrix and adjusts the genes' expression
#' for cis-eQTLs. Should be used with the 'probe_subset' feature to e.g. adjust
#' only the expression of CpG genes.
#'
#' @param expr Gene expression matrix with all the genes in the column.
#' Column names need to be illumina probe ids.
#' @param eqtls List of eqtls for which to adjust the data
#' @param geno_data Available genotype data
#' @param probe_subset Subset of probes for which to adjust the data
#'
#' @return The gene expression matrix adjusted for cis-eQTLs
#'
#' @author  Johann Hawe
#'
# ------------------------------------------------------------------------------
adjust_cis_eqtls <- function(expr, eqtls, geno_data, gene_sub=NULL){

  # sanity check
  if(!is.null(gene_sub) && !all(gene_sub %in% colnames(expr))) {
    stop("Invalid probe subset provided!")
  }

  if(!is.null(gene_sub)) {
    # for the subset of probes adjust for potential cis-eqtls (i.e. snp-effects)
    toUse <- which(eqtls$gene %in% gene_sub)
    if(length(toUse) == 0){
      return(expr)
    }
  } else {
    # use all probes
    toUse <- which(eqtls$gene %in% colnames(expr))
    if(length(toUse) == 0){
      return(expr)
    }
  }

  eqtls <- eqtls[toUse, , drop=F]
  # for each gene, check whether we have an associated cis-eqtl
  temp <- lapply(colnames(expr), function(x) {
    # indices of ciseqtls
    idxs <- which((eqtls$gene == x) & eqtls$snp_id %in% colnames(geno_data))
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
load_eqtls <- function(fceqtl, genes_sub=NULL) {
  # columns: top SNP KORA F4;Chr SNP;minor allele KORA F4;MAF KORA F4;Probe_Id;Gene;Chr Gene; etc...
  eqtls <- read.csv(fceqtl,
                    sep=";", header=T, stringsAsFactors=F)

  sids <- eqtls$top.SNP.KORA.F4
  keep <- which(rep(T,length(sids)))
  if(!is.null(genes_sub)){
    keep <- which(eqtls$Gene %in% genes_sub)
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
