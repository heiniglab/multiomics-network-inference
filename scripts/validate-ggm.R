# ------------------------------------------------------------------------------
#'
#' This script is used to validate already calculated ggm fits.
#' The fitted ggm models are loaded and validated using the concept of individual
#' link types, i.e. cpg-gene links, snp-gene links and gene-gene links.
#' @autor Johann Hawe
#'
# ------------------------------------------------------------------------------
log <- file(snakemake@log[[1]], open="wt")
sink(log)

# define easy concatenation operator
`%+%` = paste0

# ------------------------------------------------------------------------------
# Load libraries and source scripts
# ------------------------------------------------------------------------------
library(BDgraph)
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

# set some ggplot defaults
ggplot2::theme_set(theme_bw())

# ------------------------------------------------------------------------------
# Get snakemake parameters
# ------------------------------------------------------------------------------

# inputs
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
fciseqtl_joehanes <- snakemake@input[["cis_joehanes"]]
ftranseqtl_joehanes <- snakemake@input[["trans_joehanes"]]

# params
threads <- snakemake@threads
sentinel <- snakemake@wildcards$sentinel
dmediation_plots <- snakemake@params$dmediation_plots
if(!dir.exists(dmediation_plots)) {
	dir.create(dmediation_plots)
}
mediation_cutoff <- snakemake@params$mediation_cutoff

# this is the main outfile, to which to write all the validation results
fout <- snakemake@output[[1]]
fout_mediation_detail <- snakemake@output[["mediation_detail"]]

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

print("Loading gene expression data.")

# load GEO/ARCHS4  data
geo <- fread(fgeo,
             header=T, sep="\t")
colnames(geo)[1] <- "symbol"
setkey(geo, symbol)

print("Loading KORA eQTL.")
kora_ceqtl <- fread(fciseqtl_kora, sep="\t",
                    header=T)
kora_teqtl <- fread(ftranseqtl_kora, sep="\t",
                    header=T)

print("Loading Joehanes eQTL.")
# load cis-eQTL data
ceqtl <- fread(paste0("zcat ", fciseqtl_joehanes),
               sep=",")
# get only cis eqtls defined in the paper
ceqtl <- ceqtl[ceqtl$Is_Cis==1]
print(paste0("Loaded ", nrow(ceqtl), " cis-eQTL"))

# load the trans-eQTL dataset
teqtl <- fread(paste0("zcat ", ftranseqtl_joehanes),
               sep=",")

# sanity check
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
          "spath", "spath_selected", "cross_cohort_mcc", "cross_cohort_mcc_frac",
	  "mediation_total", "mediation_best_gene", "mediation_best_corr",
          "mediation_notselected_significant", "mediation_selected_significant",
          "mediation_notselected_significant.list", "mediation_selected_significant.list",
          "mediation_notselected","mediation_selected","log10_mediation",
	  "mediation_cross_cohort_correlation", "mediation_cross_cohort_fraction",
	  "mediation_cross_cohort_fraction_validation_significant",
          "cisEqtl", "transEqtl_cgenes","transEqtl_tfs", "bonder_cis_eQTM",
          "geo_gene_gene", "cohort_gene_gene",
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

  # we get the according dataset on which to validate
  # data2 is only used for the mediation analysis check,
  # where we investigate whether there would be any effect
  # in the other cohort, in case we didnt pick it up with our models
  if("lolipop" %in% cohort){
    # ggm calculated on lolipop, validate on kora
    data <- kdata
    data2 <- ldata
  } else {
    # assume kora cohort, validate on lolipop
    data <- ldata
    data2 <- kdata
  }
  # dnodes -> full set of possible nodes
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
    scluster <- paste(cl$csize, collapse=",")

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

    # retain only nodes in the largest connected component
    keep <- cl$membership==which.max(cl$csize)
    keep <- names(cl$membership[keep])

    # the nodes retained in the fitted graph model
    g <- graph::subGraph(keep, g)
    gnodes <- graph::nodes(g)
    if(!sentinel %in% gnodes) {
      g <- graph::addNode(sentinel, g)
    }

    # for now, filter only for those nodes for which data was available
    # should change once we update the data matrix...
    gnodes <- gnodes[gnodes %in% dnodes]
    # get names of all entities, total and selected by ggm
    snp <- names(ranges$sentinel_range)
    data[,snp] <- as.integer(as.character(data[,snp]))
    data2[,snp] <- as.integer(as.character(data2[,snp]))
    cpgs <- intersect(dnodes, names(ranges$cpgs))
    cpgs_selected <- cpgs[cpgs %in% gnodes]
    all_genes <- dnodes[!grepl("^rs|^cg", dnodes)]
    sgenes <- intersect(dnodes, ranges$snp_genes$SYMBOL)
    if(!is.na(snp_cluster)){
      sgenes_selected <- sgenes[sgenes %in% unlist(adj(g,snp))]
    } else {
      sgenes_selected <- c()
    }

    # cpg genes and TFs
    cgenes <- intersect(dnodes, ranges$cpg_genes$SYMBOL)
    # TODO also check TF to cpg-gene association..
    tfs <- intersect(dnodes, ranges$tfs$SYMBOL)

    # selected cpg genes and TFs
    if(length(cpgs_selected)>0) {
      cgenes_selected <- cgenes[cgenes %in% unlist(adj(g,cpgs_selected))]
      tfs_selected <- tfs[tfs %in% unlist(adj(g,cpgs_selected))]
    } else {
      cgenes_selected <- c()
      tfs_selected <- c()
    }

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

    # --------------------------------------------------------------------------
    # Check replication of complete graph across the two cohorts
    # Here we check via the MCC how well the graph structure correponds between
    # the two applications on the different cohorts
    # --------------------------------------------------------------------------
    if("lolipop" %in% cohort) {
	    cohort2 <- "kora"
    } else {
	    cohort2 <- "lolipop"
    }

    # get graph fit on other cohort
    g2 <- fits[[cohort2]][[gn]]
    # compare with largest connected component only
    ig = igraph::graph_from_graphnel(g2)
    cl = clusters(ig)
    keep <- cl$membership==which.max(cl$csize)
    keep <- names(cl$membership[keep])
    g2 <- graph::subGraph(keep, g2)
    if(!sentinel %in% graph::nodes(g2)) {
	    g2 <- graph::addNode(sentinel, g2)
    }
    # get adjacency matrices
    g_adj <- as(g, "matrix")
    g2_adj <- as(g2, "matrix")

    # ensure that we have the same nodes only in all graphs.
    # this might somewhat change results, but otherwise we
    # cant compute the MCC properly.
    use <- intersect(colnames(g_adj), colnames(g2_adj))
    g_adj <- g_adj[use,use]
    g2_adj <- g2_adj[use,use]

    # calculate performance using the DBgraph method compare()
    mcc <- compare(g_adj,g2_adj)["MCC", "estimate"]
    # the fraction of nodes retained in the overlap w.r.t. to the
    # total number of possible nodes
    mcc_frac <- ncol(g_adj)/ncol(data)
    row <- c(row, mcc, mcc_frac)

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
    med <- mediation(data, snp, sgenes, cpgs, fout_mediation_detail)
    row <- c(row, mediation.summary(med, sgenes, sgenes_selected, mediation_cutoff))

    # we also check the correspondence of the correlation values for all genes
    # in the other cohort
    med2 <- mediation(data2, snp, sgenes, cpgs)
    fout <- file.path(dmediation_plots, paste0(paste(sentinel,
						     cohort,
						     gn,
						     sep="_"),
					       ".pdf"))
    row <- c(row, compare_mediation_results(sentinel, med, med2,
			      sgenes_selected, mediation_cutoff,
			      fout))

    # --------------------------------------------------------------------------
    # (2) check cis-eQTL in independent study
    # --------------------------------------------------------------------------

    # filter ceqtl to be only related to our sentinel SNP
    # TODO use proxy/high ld snps to increase ceqtl number?
    ceqtlsub <- ceqtl[ceqtl$Rs_ID %in% snp]
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
    teqtlsub <- teqtl[teqtl$Rs_ID %in% snp]# & (teqtl$log10FDR) < (-2),,drop=F]
    if(nrow(teqtlsub)<1) {
      warning("Sentinel " %+% sentinel %+% " not available in trans-eQTL data.")
      # report NA in stats file
      row <- c(row, NA, NA)
    } else {
      row <- c(row, validate.trans.genes(teqtlsub, cgenes, tfs,
                                         cgenes_selected, tfs_selected))
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

    # load geo whole blood expression data
    geo_sym <- all_genes[all_genes %in% geo$symbol]
    geosub <- geo[geo_sym]
    geosub <- t(geosub[,-1,with=F])
    colnames(geosub) <- geo_sym

    # We collected the expression data from 3 different datasets, i.e. kora/lolipop,
    # GEO and the Geuvadis datasets.
    # Now, for each gene we calculate the correlation against each other gene in each set
    # of samples.

    # list of all expr datasets
    expr.data <- list(geo=geosub, cohort=data[,all_genes,drop=F])
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
  stop(paste0("Error for sentinel with id",
              sentinel,
              ". Not all models were validated successfully."))
}
