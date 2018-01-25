#'
#' This script is used to validate already calculated ggm fits.
#' The fitted ggm models are loaded and validated using the concept of individual
#' link types, i.e. cpg-gene links, snp-gene links and gene-gene links.
#' @autor Johann Hawe
#'

library(graph)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(qvalue)
library(Homo.sapiens)

source("R/validation.R")
source("R/lib.R")

args <- commandArgs(trailingOnly=T)
sentinel <- args[1]
if(is.null(sentinel) | sentinel=="") stop("No sentinel provided!")
#sentinel <- "rs9859077"
#sentinel <- "rs730775"

cat("Using sentinel",sentinel, "\n")

# define easy concatenation operator
`%+%` = paste0

# load preprocessed cohort data
env <- new.env()
load("results/current/data.processed.RData", envir=env)

# load external data
geu <- fread("data/current/geuvadis/GD660.GeneQuantRPKM.tsv", header=T,
             data.table=F, sep="\t")
rownames(geu) <- geu[,1]
geu <- t(geu[,-1])
geo <- fread("data/current/archs4/whole_blood/expression_matrix_norm_peer.tsv", 
             data.table=F, header=T, sep="\t")
rownames(geo) <- geo[,1]
geo <- t(geo[,-1])

kora.ceqtl <- fread("data/current/cohorts/eqtls/kora-cis-eqtls.csv", sep=";",
                         header=T, data.table=F)
kora.teqtl <- fread("data/current/cohorts/eqtls/kora-trans-eqtls.csv", sep=";",
                    header=T, data.table=F)

# load cis-eQTL data
#ceqtl <- fread("zcat data/current/joehanes-et-al-2017/eqtls/eqtl-gene-annot_cis-only.txt.gz",
#               sep=",", 
#               data.table=F, select = c("SNP_Chr","SNP_Pos_hg19","Rs_ID","SNP_MAF",
#                                        "SNP_Imputation_RSq","ProbesetID","Transcript_Chr",
#                                        "Transcript_GeneSymbol","Transcript_EntrezGeneID",
#                                        "Is_Cis","log10FDR"))
# get only cis eqtls defined in the paper
#ceqtl <- ceqtl[ceqtl$Is_Cis==1 & ceqtl$log10FDR < -2,,drop=F]
ceqtl <- kora.ceqtl
cat("Loaded ", nrow(ceqtl), " cis-eQTL\n")
head(ceqtl)

# load the trans-eQTL dataset
#teqtl <- fread("zcat data/current/joehanes-et-al-2017/eqtls/eqtl-gene-annot_trans-only_logFDR2.txt.gz", 
#               sep=",", 
#               data.table=F, select = c("SNP_Chr","SNP_Pos_hg19","Rs_ID","SNP_MAF","SNP_Imputation_RSq",
#                                        "ProbesetID","Transcript_Chr","Transcript_GeneSymbol",
#                                        "Transcript_EntrezGeneID","Is_Cis","log10FDR"))
# sanity check
#cat("only trans: ", all(teqtl$Is_Cis == 0), "\n")
teqtl <- kora.teqtl
cat("Loaded ", nrow(teqtl), " trans-eQTL\n")
head(teqtl)

# load the bonder cis-eQTMs for cpg-gene validation
# TODO we could use those for prior definitions, then use the ROADMAP 
# chromhmm states for validation...
eqtms <- fread("data/current/bonder-et-al-2017/2015_09_02_cis_eQTMsFDR0.05-CpGLevel.txt", data.table=F)

# this is the main outfile, to which to write all the validation results
ofile <- paste0("results/current/validation_nopriors/", sentinel, ".validation.txt")
# add a header line
cols <- c("sentinel","cohort","snp_genes", "snp_genes_selected", "snp_genes.list", 
          "snp_genes_selected.list", "cpg_genes", "cpg_genes_selected", "tfs", "tfs_selected",
          "spath", "spath_selected",
          "mediation_significant", "mediation_selected_significant",
          "mediation_significant.list", "mediation_selected_significant.list",
          "mediation","mediation_selected","log10_mediation",
          "cisEqtl", "transEqtl_cgenes","transEqtl_tfs", "bonder_cis_eQTM",
          "geuvadis_gene_gene", "geo_gene_gene", "cohort_gene_gene",
          "go_ids", "go_terms", "go_pvals", "go_qvals")
# the result table
tab <- cbind(t(cols))
colnames(tab) <- cols
tab <- tab[-1,]

# load both kora and lolipop data (for gene mediation)
kdata <- with(env, data[[sentinel]]$kora$data)
ldata <- with(env, data[[sentinel]]$lolipop$data)

# process both cohorts
cohorts <- c("lolipop", "kora")

temp <- lapply(cohorts, function(cohort){
  
  # write basic info to output
  row <- c(sentinel, cohort)
  
  ## Data preparation
  # Before starting prepare the needed data, i.e. the ggm fit, the ranges originally
  # used and the original data.matrix. We als retrieve the different set of entities,
  # i.e. the snp, the cpgs as well as the respective locus genes, tfs and shortest path
  # genes. Those entities can either have been selected via ggm graph or not.
  
  load(paste0("results/current/fits_nopriors/", sentinel, ".", cohort, ".RData"))
  
  # dnodes -> full set of possible nodes
  if("lolipop" %in% cohort){
    # ggm calculated on lolipop, validate on kora
    data <- kdata
  } else {
    # assume kora cohort, validate on lolipop
    data <- ldata
  }
  dnodes <- colnames(data)
  
  # the nodes retained in the fitted graph model
  gnodes <- nodes(graph)
  
  # for now, filter only for those nodes for which data was available
  # should change once we update the data matrix...
  gnodes <- gnodes[gnodes %in% dnodes]
  # get names of all entities, total and selected by ggm
  snp <- names(ranges$sentinel.range)
  data[,snp] <- as.integer(as.character(data[,snp]))
  cpgs <- intersect(dnodes, names(ranges$cpgs))
  cpgs.selected <- cpgs[cpgs %in% gnodes]
  all.genes <- colnames(data)[!grepl("^rs|^cg", colnames(data))]
  sgenes <- intersect(dnodes, ranges$snp.genes$SYMBOL)
  sgenes.selected <- sgenes[sgenes %in% unlist(adj(graph,snp))]
  cgenes <- intersect(dnodes, ranges$cpg.genes$SYMBOL)
  cgenes.selected <- cgenes[cgenes %in% unlist(adj(graph,cpgs))]
  
  # TODO also check TF to cpg-gene association..
  tfs <- intersect(dnodes, ranges$tfs$SYMBOL)
  tfs.selected <- tfs[tfs %in% unlist(adj(graph,cpgs))]
  tfs.selected
  
  # the shortest path genes
  spath <- ranges$spath$SYMBOL
  spath.selected <- spath[spath %in% gnodes]
  
  # write to stats file
  row <- c(row,
           length(sgenes), length(sgenes.selected),
           paste(sgenes, collapse=";"), paste(sgenes.selected, collapse=";"),
           length(cgenes), length(cgenes.selected),
           length(tfs), length(tfs.selected), 
           length(spath), length(spath.selected))
  
  # write to user
  cat("Summary on number of genes, total vs selected via ggm:\n")
  cat("sgenes\t", length(sgenes),"\t", length(sgenes.selected), "\n")
  cat("cgenes\t", length(cgenes),"\t", length(cgenes.selected), "\n")
  cat("tfs\t", length(tfs),"\t", length(tfs.selected), "\n")
  cat("spath\t", length(spath),"\t", length(spath.selected), "\n")
  
  # Above the number of snp genes (sgenes), cpg genes (cgenes), transcription
  # factors (tfs) as well as the genes on the shortest paths between tfs and sgenes (spath)
  # are shown. Starting from those set definitions, we now perform the validation of the model.
  # 
  ## SNP-gene validation
  # For the SNP-gene validation for which we have several approaches:
  # 
  # 1. Perform gene mediation. 
  # Here we check for each gene individually, whether it likely
  # mediates the effect the sentinel snp has on the defined trans-cpgs. We calculate linear models
  # between snp/genes, snp/cpgs and genes/cpgs. The $\beta$ of the linear models should add up, i.e.
  # $\beta_{sc} = \beta_{sg}*\beta_{gc}$. For each gene we get as a result the correlation between
  # the $\beta_{sc}$ and $\beta_{sg}*\beta_{gc}$ over all CpGs.
  # 2. Check cis-eQTLs. Here we check whether the selected SNP genes (connected to the sentinel) 
  # are cis-eQTLs in an independent study (Joehanes et al. 2017). The more selected genes are cis-eQTLs
  # as compared to the not selected ones, the more significant is the selection of genes in our models.
  # 3. check whether SNPs are in TFBS (in case of >1 SNPs)
  # 4. ChromHMM states for SNPs (in case of >1 SNPs)
  # concept.
  # 
  
  # (1) Perform mediation analysis
  # mediation over all snp genes
  med <- mediation(data, snp, sgenes, cpgs)
  # mediation for only ggm selected snp genes
  med.selected <- med[sgenes.selected]
  med.selected.sign <- med.selected[unlist(lapply(med.selected, 
                                                  function(x) { 
                                                    x["pvalue"]<0.05 
                                                  }))]
  
  # mediation for only not selected snp genes
  med <- med[setdiff(sgenes,sgenes.selected)]
  med.sign <- med[unlist(lapply(med, 
                                function(x) { 
                                  x["pvalue"]<0.05 
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
  cat("Difference (log10(ns/s)): ", log10(med.pv/med.pv.selected), "\n")
  # add to row
  row <- c(row,
           length(med.sign),
           length(med.selected.sign),
           paste(names(med.sign), collapse=";"),
           paste(names(med.selected.sign), collapse=";"),
           med.pv, 
           med.pv.selected,
           log10(med.pv/med.pv.selected))
  
  # (2) check cis-eQTL in independent study
  
  # filter ceqtl to be only related to our sentinel SNP
  # TODO use proxy/high ld snps to increase ceqtl number?
  ceqtlsub <- ceqtl[ceqtl$Rs_ID %in% snp,,drop=F]
  if(nrow(ceqtlsub) < 1) {
    warning("Sentinel " %+% sentinel %+% " not found in cis-eQTL data")
    # report NA in stats file 
    row <- c(row, NA)
  } else {
    ceqtl.sgenes <- sgenes[sgenes %in% ceqtlsub$Transcript_GeneSymbol]
    ceqtl.sgenes.selected <- intersect(ceqtl.sgenes, sgenes.selected)
    ceqtl.sgenes
    ceqtl.sgenes.selected
  
    # create matrix for fisher test
    cont <- matrix(c(length(ceqtl.sgenes),length(ceqtl.sgenes.selected),
                     length(sgenes),length(sgenes.selected)),
                   nrow=2,ncol=2, byrow = T)
    rownames(cont) <- c("ceqtl", "no ceqtl")
    colnames(cont) <- c("not selected", "selected")
    cont
    
    row <- c(row,
             fisher.test(cont)$p.value)
  }
  
  # (1) load chromHMM annotation (SNP annotation)
  ## not yet needed
  
  
  # (3) SNPs in TFBS
  ## not yet needed
  
  ## CpG-gene and TF-CpG validation
  # The next step is validating the cpg-gene and cpg-tf links. For this we use 3 approaches:
  # 1. check trans-eQTL in independent dataset. Here we found the paper of [Joehanes et al. 2017
  # in Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1142-6)
  # 2. check epigenetic annotation (implies change in currently performed prior calculation)
  # 3. check CpGs in TFBS (split cell-lines for priors/validation).
  
  # filter teqtl to be only related to our sentinel SNP
  # TODO use proxy/high ld snps to increase teqtl number?
  teqtlsub <- teqtl[teqtl$`top SNP KORA` %in% snp,,drop=F]# & (teqtl$log10FDR) < (-2),,drop=F]
  if(nrow(teqtlsub)<1) {
    warning("Sentinel " %+% sentinel %+% " not available in trans-eQTL data.")
    # report NA in stats file 
    row <- c(row,
             NA,
             NA)
  } else {
    teqtl.cgenes <- cgenes[cgenes %in% teqtlsub$Transcript_GeneSymbol]
    teqtl.cgenes.selected <- intersect(teqtl.cgenes, cgenes.selected)
  
    cat("CpG genes: ", cgenes, "\n")
    cat("CpG genes with trans-eQTL: ", teqtl.cgenes, "\n")
    cat("selected CpG genes with trans-eQTL: ", teqtl.cgenes.selected, "\n")
  
    # create matrix for fisher test
    cont <- matrix(c(length(teqtl.cgenes),length(teqtl.cgenes.selected),
                     length(cgenes),length(cgenes.selected)),
                   nrow=2,ncol=2, byrow = T)
    cat("confusion matrix for cgenes:\n")
    rownames(cont) <- c("teqtl", "no teqtl")
    colnames(cont) <- c("not selected", "selected")
    cont
    # report cpg-gene fisher test
    row <- c(row,
             fisher.test(cont)$p.value)
    
    # analyze the tfs, total and selected
    teqtl.tfs <- tfs[tfs %in% teqtlsub$Transcript_GeneSymbol]
    teqtl.tfs.selected <- intersect(teqtl.tfs, tfs.selected)
  
    cat("TFs: ", tfs, "\n")
    cat("TFs with trans-eQTL: ", teqtl.tfs, "\n")
    cat("selected TFs with trans-eQTL: ", teqtl.tfs.selected, "\n")
    # create matrix for fisher test
    cont <- matrix(c(length(teqtl.tfs),length(teqtl.tfs.selected),
                     length(tfs),length(tfs.selected)),
                   nrow=2,ncol=2, byrow = T)
    cat("confusion matrix for TFs:\n")
    rownames(cont) <- c("teqtl", "no teqtl")
    colnames(cont) <- c("not selected", "selected")
    cont
    # report tf fisher test
    row <- c(row,
             fisher.test(cont)$p.value)
  }
  
  # we can also count how many of the identified cpg
  # genes are eQTMs in the bonder dataset
  
  # first convert bonder ensembl ids to symbols
  # (bonder results already filtered for sign. assoc.)
  bnd.symbol <- symbol.from.ensembl(eqtms$ProbeName)$SYMBOL
  bnd.symbol <- bnd.symbol[!is.na(bnd.symbol)]
  v1 <- sum(setdiff(cgenes, cgenes.selected) %in% bnd.symbol)
  v2 <- sum(!setdiff(cgenes, cgenes.selected) %in% bnd.symbol)
  v3 <- sum(cgenes.selected %in% bnd.symbol)
  v4 <- sum(!cgenes.selected %in% bnd.symbol)
  # build the confusion matrix 
  cat("confusion matrix for cpg-genes:\n")
  cont <- matrix(c(v3,v1,v4,v2), 
                 ncol=2, byrow=T)
  rownames(cont) <- c("has_eQTM", "has_no_eQTM")
  colnames(cont) <- c("ggm_selected", "not_gmm_selected")
  cont
  # report fisher test
  row <- c(row,
           fisher.test(cont)$p.value)

  ## Gene-Gene validation
  # Finally we perform the gene-gene validation:
  # 1. check co-expression in independent data. 
  # 2. gene set enrichment
  # 3. check co-citation 
  # 
  # For the co-expression analysis we load the data from Geuvadis (Lappaleinen et al.),
  # this is gene expression data from 462 LCL cell-lines. The data already is normalized
  # to RPKMs and technical artifacts are removed (using PEER factors). We currently
  # calculated correlations only on the subset of genes which is within our network.
  # load the expression data from geuvadis/lappaleinen
 
  # subset to genes which were available to the gggm algorithm
  geusub <- geu[,colnames(geu) %in% dnodes, drop=F]
  
  # load geo whole blood expression data
  geosub <- geo[,colnames(geo) %in% dnodes, drop=F]
  
  # We collected the expression data from 3 different datasets, i.e. kora/lolipop,
  # GEO and the Geuvadis datasets.
  # Now, for each gene we calculate the correlation against each other gene in each set
  # of samples. 
  
  # get the gene-gene links which are possible in our graph
  em <- expand.grid(all.genes, all.genes, stringsAsFactors = F)
  em <- (em[em$Var1 != em$Var2, ])
  # remove duplicated entries
  em.names <- vector(mode="character", length=nrow(em))
  for(i in 1:nrow(em)){
    em.names[i] <- paste0(sort(c(em[i,1], em[i,2])), collapse="")
  }
  em <- em[!duplicated(em.names),]
  em <- cbind(em[,1], em[,2])
  em <- as.data.frame(em, stringsAsFactors=F)
  em
  
  # get the gene-gene links which can be found in our graph
  em.selected <- t(edgeMatrix(graph))
  em.selected <- cbind(gnodes[em.selected[,1]], gnodes[em.selected[,2]])
  em.selected <- em.selected[!grepl("^rs|^cg", em.selected[,1]),]
  em.selected <- em.selected[!grepl("^rs|^cg", em.selected[,2]),]
  colnames(em) <- colnames(em.selected) <- c("n1", "n2")
  em.selected <- as.data.frame(em.selected, stringsAsFactors=F)
  em.selected
  
  # list of all expr datasets
  expr.data <- list(geuvadis=geusub, geo=geosub, cohort=data[,!grepl("^cg|^rs",colnames(data)),drop=F])
  
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
  
  # report fisher test result for each DS
  row <- c(row,
           results[[1]],
           results[[2]],
           results[[3]])
  
  ## GO enrichment
  # do one for all gene nodes in the graph
  source("R/go-enrichment.R")
  gn <- gnodes[!grepl("^rs|^cg",gnodes)]
  
  # define background set
  # for now all possible symbols from the array annotation
  # are used
  #bgset <- dnodes[!grepl("^rs|^cg",dnodes)]
  require(illuminaHumanv3.db)
  annot <- illuminaHumanv3SYMBOLREANNOTATED
  bgset <- unique(unlist(as.list(annot)))
  
  if(length(gn)<length(bgset)){
    go.tab <- go.enrichment(gn, bgset, gsc)
    go.tab <- go.tab[go.tab$q<0.01,,drop=F]
    if(nrow(go.tab)>0){
      row <- c(row,
               paste0(go.tab$GOID, collapse=","),
               paste0(go.tab$Term, collapse=","),
               paste0(go.tab$Pvalue, collapse=","),
               paste0(go.tab$q, collapse=","))
    } else {
      row <- c(row,
               NA,
               NA,
               NA,
               NA)
    }
    
  }
  # finish the current sentinel
  # all went well
  if(length(row) == ncol(tab)){
    tab <<- rbind(tab,row)
  }
})

# report results
if(nrow(tab) == 2){
  # write output file
  write.table(file=ofile, tab, col.names=T, row.names=F, quote=F,sep="\t")
} else {
  # report error
  cat(file=ofile, "Error for sentinel with id", sentinel, "\n")
}

