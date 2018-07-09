#'
#' This script is used to validate already calculated ggm fits.
#' The fitted ggm models are loaded and validated using the concept of individual
#' link types, i.e. cpg-gene links, snp-gene links and gene-gene links.
#' @autor Johann Hawe
#'

# define easy concatenation operator
`%+%` = paste0

# ------------------------------------------------------------------------------
# Load libraries and source scripts
# ------------------------------------------------------------------------------
library(igraph)
library(graph)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(qvalue)
library(Homo.sapiens)
library(illuminaHumanv3.db)

source("scripts/go-enrichment.R")
source("scripts/validation.R")
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# Get snakemake parameters
# ------------------------------------------------------------------------------

sentinel <- snakemake@wildcards$sentinel
fkora_data <- snakemake@input[["kora_data"]]
franges <- snakemake@input[["ranges"]]
flolipop_data <- snakemake@input[["lolipop_data"]]
fkora_fit <- snakemake@input[["kora_fit"]]
flolipop_fit <- snakemake@input[["lolipop_fit"]]
fgtex <- snakemake@input[["gtex"]]
fgeo <- snakemake@input[["geo"]]
fciseqtl_kora <- snakemake@input[["cis_kora"]]
ftranseqtl_kora <- snakemake@input[["trans_kora"]]
fbonder_eqtm <- snakemake@input[["bonder_eqtm"]]

# this is the main outfile, to which to write all the validation results
fout <- snakemake@output[[1]]

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

print("Loading gene expression data.")

# load gtex data
geu <- fread(fgtex, header=T,
             data.table=F, sep="\t")
rownames(geu) <- geu[,1]
geu <- t(geu[,-1])
geo <- fread(fgeo, 
             data.table=F, header=T, sep="\t")
rownames(geo) <- geo[,1]
geo <- t(geo[,-1])

print("Loading KORA eQTL.")

kora_ceqtl <- fread(fciseqtl_kora, sep=";",
                         header=T, data.table=F)
kora_teqtl <- fread(ftranseqtl_kora, sep=";",
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
ceqtl <- kora_ceqtl
print(paste0("Loaded ", nrow(ceqtl), " cis-eQTL"))

# load the trans-eQTL dataset
#teqtl <- fread("zcat data/current/joehanes-et-al-2017/eqtls/eqtl-gene-annot_trans-only_logFDR2.txt.gz", 
#               sep=",", 
#               data.table=F, select = c("SNP_Chr","SNP_Pos_hg19","Rs_ID","SNP_MAF","SNP_Imputation_RSq",
#                                        "ProbesetID","Transcript_Chr","Transcript_GeneSymbol",
#                                        "Transcript_EntrezGeneID","Is_Cis","log10FDR"))
# sanity check
#cat("only trans: ", all(teqtl$Is_Cis == 0), "\n")
teqtl <- kora_teqtl
print(paste0("Loaded ", nrow(teqtl), " trans-eQTL"))

# load the bonder cis-eQTMs for cpg-gene validation
# TODO we could use those for prior definitions, then use the ROADMAP 
# chromhmm states for validation...
eqtms <- fread(fbonder_eqtm, data.table=F)

# add a header line
cols <- c("sentinel","cohort","graph_type", "number_nodes", "number_edges", 
          "graph_density", "cluster", "cluster_sizes", 
          "snp_cluster", "snp_genes", "snp_genes_selected", "snp_genes.list", 
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

print("Loading cohort data.")

# load both kora and lolipop data (for gene mediation)
kdata <- readRDS(fkora_data)
ldata <- readRDS(flolipop_data)

ranges <- readRDS(franges)

# rds files containing ggm fits
print("Loading GGM fits.")

kfit <- readRDS(fkora_fit)
lfit <- readRDS(flolipop_fit)
fits <- list(kora=kfit, lolipop=lfit)

# process both cohorts
cohorts <- c("lolipop", "kora")

# ------------------------------------------------------------------------------
# Main apply for validation, apply for both cohorts
# ------------------------------------------------------------------------------
temp <- lapply(cohorts, function(cohort){
  
  print(paste0("Validating ", cohort, " fits. -------------------------------"))
        
  # write basic info to output
  row <- c(sentinel, cohort)

  # ----------------------------------------------------------------------------  
  ## Data preparation
  # ----------------------------------------------------------------------------
  # Before starting prepare the needed data, i.e. the ggm fit, the ranges originally
  # used and the original data.matrix. We als retrieve the different set of entities,
  # i.e. the snp, the cpgs as well as the respective locus genes, tfs and shortest path
  # genes. Those entities can either have been selected via ggm graph or not.
  graph <- fits[[cohort]]$graph
  graph_no_priors <- fits[[cohort]]$graph_no_priors
  graph_genenet <- fits[[cohort]]$graph_genenet
  
  # we apply over all graphs...
  gs <- listN(graph, graph_no_priors, graph_genenet)
  
  # dnodes -> full set of possible nodes
  if("lolipop" %in% cohort){
    # ggm calculated on lolipop, validate on kora
    data <- kdata
  } else {
    # assume kora cohort, validate on lolipop
    data <- ldata
  }
  dnodes <- colnames(data)
  
  # ----------------------------------------------------------------------------
  # we have now all we need, now get specific for the 
  # individual graphs (fittet with and w/o priors, genenet...)
  # ----------------------------------------------------------------------------
  rows <- lapply(names(gs), function(gn) {

    print(paste0("Validating ", gn, " graph fit. ----------------------------"))
    
    # --------------------------------------------------------------------------
    # remember whether currently we are looking at the
    # graph with or without priros by using the name
    # --------------------------------------------------------------------------
    row <- c(row, gn)

    # get the graph
    g <- gs[[gn]]
    
    # get basic stats (number nodes, number edges, graph density)
    nn <- numNodes(g)
    ne <- numEdges(g)
    gd <- (ne * 2) / (nn * (nn-1))
    
    # --------------------------------------------------------------------------
    # add basic stats to row
    # --------------------------------------------------------------------------
    row <- c(row, nn, ne, gd)
    
    # check clusters and get the cluster with the SNP
    ig = igraph::graph_from_graphnel(g)
    cl = clusters(ig)
    ncluster <- cl$no
    scluster <- paste(sort(cl$csize, decreasing = T), 
                      collapse=",")
    
    # remember snp membership
    snp_cluster <- NA
    if(sentinel %in% names(cl$membership)) {
      snp_cluster <- cl$membership[sentinel]
      snp_cluster_size <- cl$csize[snp_cluster]
      
      snp_cluster <- paste0(snp_cluster, "-", snp_cluster_size)
    } 
    # --------------------------------------------------------------------------
    # add cluster information to row
    # --------------------------------------------------------------------------
    row <- c(row, ncluster, scluster, snp_cluster)
    
    # the nodes retained in the fitted graph model
    gnodes <- nodes(g)
    
    # for now, filter only for those nodes for which data was available
    # should change once we update the data matrix...
    gnodes <- gnodes[gnodes %in% dnodes]
    # get names of all entities, total and selected by ggm
    snp <- names(ranges$sentinel_range)
    data[,snp] <- as.integer(as.character(data[,snp]))
    cpgs <- intersect(dnodes, names(ranges$cpgs))
    cpgs_selected <- cpgs[cpgs %in% gnodes]
    all_genes <- colnames(data)[!grepl("^rs|^cg", colnames(data))]
    sgenes <- intersect(dnodes, ranges$snp_genes$SYMBOL)
    if(!is.na(snp_cluster)){
      sgenes_selected <- sgenes[sgenes %in% unlist(adj(g,snp))]
    } else {
      sgenes_selected <- c()
    }
    
    cgenes <- intersect(dnodes, ranges$cpg_genes$SYMBOL)
    cgenes_selected <- cgenes[cgenes %in% unlist(adj(g,cpgs_selected))]
    
    # TODO also check TF to cpg-gene association..
    tfs <- intersect(dnodes, ranges$tfs$SYMBOL)
    tfs_selected <- tfs[tfs %in% unlist(adj(g,cpgs_selected))]
    tfs_selected
    
    # the shortest path genes
    spath <- ranges$spath$SYMBOL
    spath_selected <- spath[spath %in% gnodes]
    
    # --------------------------------------------------------------------------
    # write to stats file
    # --------------------------------------------------------------------------
    row <- c(row,
             length(sgenes), length(sgenes_selected),
             paste(sgenes, collapse=";"), paste(sgenes_selected, collapse=";"),
             length(cgenes), length(cgenes_selected),
             length(tfs), length(tfs_selected), 
             length(spath), length(spath_selected))
    
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
    
    # --------------------------------------------------------------------------
    # (1) Perform mediation analysis
    # mediation over all snp genes
    # --------------------------------------------------------------------------
    med <- mediation(data, snp, sgenes, cpgs)
    row <- c(row, mediation.summary(med, sgenes, sgenes_selected))
    
    # --------------------------------------------------------------------------
    # (2) check cis-eQTL in independent study
    # --------------------------------------------------------------------------
    
    # filter ceqtl to be only related to our sentinel SNP
    # TODO use proxy/high ld snps to increase ceqtl number?
    ceqtlsub <- ceqtl[ceqtl$Rs_ID %in% snp,,drop=F]
    if(nrow(ceqtlsub) < 1) {
      warning("Sentinel " %+% sentinel %+% " not found in cis-eQTL data")
      # report NA in stats file 
      row <- c(row, NA)
    } else {
      ceqtl_sgenes <- sgenes[sgenes %in% ceqtlsub$Transcript_GeneSymbol]
      ceqtl_sgenes_selected <- intersect(ceqtl_sgenes, sgenes_selected)
      
      # create matrix for fisher test
      cont <- matrix(c(length(ceqtl_sgenes),length(ceqtl_sgenes_selected),
                       length(sgenes),length(sgenes_selected)),
                     nrow=2,ncol=2, byrow = T)
      rownames(cont) <- c("ceqtl", "no ceqtl")
      colnames(cont) <- c("not selected", "selected")
      cont
      # ------------------------------------------------------------------------
      # add cis eqtl summary to stats
      # ------------------------------------------------------------------------
      row <- c(row,
               fisher.test(cont)$p.value)
    }
    
    # (1) load chromHMM annotation (SNP annotation)
    ## not yet needed
    
    
    # (3) SNPs in TFBS
    ## not yet needed
    
    # --------------------------------------------------------------------------
    ## CpG-gene and TF-CpG validation
    # --------------------------------------------------------------------------
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
      row <- c(row, NA, NA)
    } else {
      row <- c(row, validate.trans.genes(teqtl, cgenes, tfs, 
                                         cgenes_selected,   tfs_selected))
    }
    
    # we can also count how many of the identified cpg
    # genes are eQTMs in the bonder dataset
    # first convert bonder ensembl ids to symbols
    # (bonder results already filtered for sign. assoc.)
    bnd_symbol <- symbol.from.ensembl(eqtms$ProbeName)$SYMBOL
    bnd_symbol <- bnd_symbol[!is.na(bnd_symbol)]
    row <- c(row, validate.cpggenes(bnd_symbol, cgenes, cgenes_selected))
  
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
    
    # list of all expr datasets
    expr.data <- list(geuvadis=geusub, geo=geosub, cohort=data[,!grepl("^cg|^rs",colnames(data)),drop=F])
    row <- c(row, validate.gene2gene(expr.data, g, all_genes))
    
    # --------------------------------------------------------------------------
    ## GO enrichment
    # do one for all gene nodes in the graph
    # --------------------------------------------------------------------------
    row <- c(row, validate.geneenrichment(gnodes))
    # finish the current sentinel
    # all went well
    if(length(row) == ncol(tab)){
      tab <<- rbind(tab,row)
    }
  })
})

# report results
if(nrow(tab) == 6){
  # write output file
  write.table(file=fout, tab, col.names=T, 
              row.names=F, quote=F, sep="\t")
} else {
  # report error
  cat(file=fout, "Error for sentinel with id", sentinel, "\n")
}
