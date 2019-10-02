# ------------------------------------------------------------------------------
#'
#' This script is used to validate already calculated ggm fits.
#' The fitted ggm models are loaded and validated using the concept of individual
#' link types, i.e. cpg-gene links, snp-gene links and gene-gene links.
#'
#'  @autor Johann Hawe
#'
# ------------------------------------------------------------------------------
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# define easy concatenation operator
`%+%` = paste0

# ------------------------------------------------------------------------------
print("Load libraries and source scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(BDgraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(illuminaHumanv3.db))
suppressPackageStartupMessages(library(ggplot2))
library(parallel)
library(reshape2)
library(cowplot)

source("scripts/go-enrichment.R")
source("scripts/validation_methods.R")
source("scripts/mediation_methods.R")
source("scripts/lib.R")
source("scripts/reg_net.R")

# ------------------------------------------------------------------------------
print("Get snakemake parameters.")
# ------------------------------------------------------------------------------

# inputs
fkora_data <- snakemake@input[["kora_data"]]
franges <- snakemake@input[["ranges"]]
flolipop_data <- snakemake@input[["lolipop_data"]]
fkora_fit <- snakemake@input[["kora_fit"]]
flolipop_fit <- snakemake@input[["lolipop_fit"]]

# contain the 'old' fits, i.e. old glasso and w/o genie3
#fkora_fit_old <- snakemake@input[["kora_fit_old"]]
#flolipop_fit_old <- snakemake@input[["lolipop_fit_old"]]

fgtex <- snakemake@input[["gtex"]]
fgeo <- snakemake@input[["geo"]]
fciseqtl_kora <- snakemake@input[["cis_kora"]]
ftranseqtl_kora <- snakemake@input[["trans_kora"]]
fbonder_eqtm <- snakemake@input[["bonder_eqtm"]]
fciseqtl_joehanes <- snakemake@input[["cis_joehanes"]]
ftranseqtl_joehanes <- snakemake@input[["trans_joehanes"]]
fmediation_kora <- snakemake@input[["mediation_kora"]]
fmediation_lolipop <- snakemake@input[["mediation_lolipop"]]

# params
threads <- snakemake@threads
sentinel <- snakemake@wildcards$sentinel
mediation_cutoff <- snakemake@params$mediation_cutoff
cohort <- snakemake@wildcards$cohort

# output
# main outfile for validation results
fout <- snakemake@output[[1]]

# additional outputs, mostly plots
#fmediation_summary_plot <- snakemake@output$mediation_summary_plot

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------

print("Loading gene expression data.")

# load GEO/ARCHS4  data
geo <- fread(fgeo,
             header = T, sep = "\t")
colnames(geo)[1] <- "symbol"
setkey(geo, symbol)

print("Loading Joehanes eQTL.")
ceqtl <- fread(fciseqtl_joehanes,
               sep = ",")
# get only cis eqtls defined in the paper
ceqtl <- ceqtl[ceqtl$Is_Cis == 1]
print(paste0("Loaded ", nrow(ceqtl), " cis-eQTL"))

# trans eQTL
teqtl <- fread(ftranseqtl_joehanes,
               sep = ",")
print(paste0("Loaded ", nrow(teqtl), " trans-eQTL"))

# load the bonder cis-eQTMs for cpg-gene validation
eqtms <- fread(fbonder_eqtm, data.table = F)

# ------------------------------------------------------------------------------
print("Loading cohort data.")
# ------------------------------------------------------------------------------
kdata <- readRDS(fkora_data)
ldata <- readRDS(flolipop_data)

# remove (rare) all-NA cases. This can happen due to scaling of all-zero entities,
# which can arise due to a very large number of cis-meQTLs which effects get
# removed from the CpGs during data preprocessing.
# NOTE: we could possibly handle this differently? Seems that these particular
# cpgs are highly influenced by genetic effects?
use <- apply(kdata,2,function(x) (sum(is.na(x)) / length(x)) < 1)
kdata <- kdata[,use]
use <- apply(ldata,2,function(x) (sum(is.na(x)) / length(x)) < 1)
ldata <- ldata[,use]

# load ranges
ranges <- readRDS(franges)

# ------------------------------------------------------------------------------
print("Loading GGM fits.")
# ------------------------------------------------------------------------------
kfit <- readRDS(fkora_fit)
lfit <- readRDS(flolipop_fit)

#kfit_old <- readRDS(fkora_fit_old)
#lfit_old <- readRDS(flolipop_fit_old)

fits <- list(kora = kfit, lolipop = lfit,
             kora_old=kfit_old, lolipop_old=lfit_old)

# we get the according dataset on which to validate
if ("lolipop" %in% cohort) {
  # ggm calculated on lolipop, validate on kora
  data_val <- kdata
  data_fit <- ldata
  cohort_val <- "kora"
} else {
  # assume kora cohort, validate on lolipop
  data_val <- ldata
  data_fit <- kdata
  cohort_val <- "lolipop"
}

graph_types <- c("bdgraph", "bdgraph_no_priors",
                 "genenet", "irafnet",
                 "glasso", "glasso_no_priors", "genie3")

# validate all graph models
# NOTE: although we used multi-threading before, it seems that this results in
# problems when integrating the GO enrichment (cryptic database disk image errors)
# Therefore we only use one thread for now, but run a lot of distinct jobs, this
# seems to do the trick
valid <- mclapply(graph_types, function(graph_type) {
  print(paste0("Validating ", cohort, " fit for '", graph_type , "' graph fit."))

  row <- c(sentinel = sentinel,
           cohort = cohort,
           graph_type = graph_type)

  # ----------------------------------------------------------------------------
  print("Preparing fit.")
  # ----------------------------------------------------------------------------
  # those were adjusted and newly fitted
 # if(graph_type %in% c("glasso", "glasso_no_priors", "genie3")) {
    graph <- fits[[cohort]][[graph_type]]
 # } else {
 #   graph <- fits[[paste0(cohort, "_old")]][[graph_type]]
 # }

  # dnodes -> full set of possible nodes
  dnodes <- colnames(data_val)

  # ----------------------------------------------------------------------------
  print("Getting basic stats (number nodes, number edges, graph density)")
  # ----------------------------------------------------------------------------
  nn <- numNodes(graph)
  ne <- numEdges(graph)
  gd <- (ne * 2) / (nn * (nn - 1))
  row <- c(row, number_nodes=nn, number_edges=ne, graph_density=gd)

  # ------------------------------------------------------------------------------
  print("Getting cluster information")
  # ------------------------------------------------------------------------------
  # get all clusters in the graph
  ig = igraph::graph_from_graphnel(graph)
  cl = clusters(ig)
  ncluster <- cl$no
  scluster <- paste(cl$csize, collapse = ",")

  # remember snp membership
  snp_cluster <- NA
  snp_cluster_size <- NA
  if (sentinel %in% names(cl$membership)) {
    snp_cluster <- cl$membership[sentinel]
    snp_cluster_size <- cl$csize[snp_cluster]
  }

  row <- c(row,
           cluster=ncluster,
           cluster_sizes=scluster,
           snp_cluster=unname(snp_cluster),
           snp_cluster_size=snp_cluster_size)

  # ----------------------------------------------------------------------------
  print("Getting largest CC for validation.")
  # ----------------------------------------------------------------------------
  keep <- cl$membership == which.max(cl$csize)
  keep <- names(cl$membership[keep])
  if(!is.null(keep)) {
    graph_maxcluster <- graph::subGraph(keep, graph)
  }

  # the nodes retained in the fitted graph model in the largest CC
  gnodes <- graph::nodes(graph_maxcluster)

  # ----------------------------------------------------------------------------
  print("Calculating graph score.")
  # ----------------------------------------------------------------------------
  # we use the (full) igraph object for this, will be filtered for the sentinel
  # cluster
  score <- get_graph_score(ig, sentinel, ranges, gd)
  row <- c(row,
           graph_score = score)

  # ----------------------------------------------------------------------------
  print("Defining entity sets (selected / not selected)")
  # ----------------------------------------------------------------------------

  # get names of all entities, total and selected by ggm
  snp <- sentinel
  data_val[, snp] <- as.integer(as.character(data_val[, snp]))
  data_fit[, snp] <- as.integer(as.character(data_fit[, snp]))

  if(ranges$seed == "meqtl") {
    trans_entities <- intersect(dnodes, names(ranges$cpgs))
  } else {
    trans_entities <- intersect(dnodes, ranges$trans_genes$SYMBOL)
  }
  trans_entities_selected <- trans_entities[trans_entities %in% gnodes]

  all_genes <- dnodes[!grepl("^rs|^cg", dnodes)]
  sgenes <- intersect(dnodes, ranges$snp_genes$SYMBOL)
  if (snp %in% gnodes) {
    sgenes_selected <- sgenes[sgenes %in% unlist(adj(graph_maxcluster, snp))]

    # only proceed if we have at least one selected snp gene and trans entity
    if(length(sgenes_selected) > 0 & length(trans_entities_selected) > 0) {

      # collection of 'real' selected SNP genes (with path to one of the trans
      # entities)
      temp_genes <- c()

      # check each snp gene individually
      for(sg in sgenes_selected) {
        # snp genes have to be connected to trans entities without traversing
        # the snp itself and must not go via another snp gene
        # -> create subgraph without those entities
        sg_other <- setdiff(sgenes_selected, sg)
        graph_temp <- subGraph(setdiff(gnodes, c(snp, sg_other)),
                               graph_maxcluster)
        igraph_temp <- igraph::graph_from_graphnel(graph_temp)

        # get paths
        paths <-
          suppressWarnings(get.shortest.paths(igraph_temp, sg,
                                              intersect(V(igraph_temp)$name,
                                                        trans_entities_selected)))$vpath[[1]]

        # did we find a path with the current SNP gene?
        temp_genes <- unique(c(temp_genes, intersect(sg, paths$name)))
      }
      sgenes_selected <- temp_genes
    } else {
      sgenes_selected <- c()
    }
  } else {
    sgenes_selected <- c()
  }

  # cpg genes and TFs
  cgenes <- intersect(dnodes, ranges$cpg_genes$SYMBOL)
  # TODO also check TF to cpg-gene association..
  tfs <- intersect(dnodes, ranges$tfs$SYMBOL)

  if(ranges$seed == "meqtl") {
    # selected cpg genes and TFs
    if (length(trans_entities_selected) > 0) {
      cgenes_selected <- cgenes[cgenes %in% unlist(adj(graph_maxcluster,
                                                       trans_entities_selected))]
    } else {
      cgenes_selected <- c()
    }
  } else {
    cgenes_selected <- c()
  }

  if(length(trans_entities_selected) > 0) {
    tfs_selected <- tfs[tfs %in% unlist(adj(graph_maxcluster,
                                            trans_entities_selected))]
  } else {
    tfs_selected <- c()
  }

  # the shortest path genes
  spath <- ranges$spath$SYMBOL
  spath_selected <- spath[spath %in% gnodes]

  # add to row
  row <- c(
    row,
    snp_genes=length(sgenes),
    snp_genes_selected=length(sgenes_selected),
    snp_genes.list=paste(sgenes, collapse = ";"),
    snp_genes_selected.list=paste(sgenes_selected, collapse = ";"),
    trans_entities = length(trans_entities),
    trans_entities_selected = length(trans_entities_selected),
    cpg_genes=length(cgenes),
    cpg_genes_selected=length(cgenes_selected),
    tfs=length(tfs),
    tfs_selected=length(tfs_selected),
    spath=length(spath),
    spath_selected=length(spath_selected),
    tfs_per_trans=length(tfs)/length(trans_entities),
    tfs_per_trans_selected=length(tfs_selected)/length(trans_entities_selected)
  )

  # ------------------------------------------------------------------------------
  print("Using MCC to check how well graph replicated across cohorts.")
  # ------------------------------------------------------------------------------
  # get graph fit on other cohort
  # those were adjusted and newly fitted
  #if(graph_type %in% c("glasso", "glasso_no_priors", "genie3")) {
    graph_val <- fits[[cohort_val]][[graph_type]]
  #} else {
  #  graph_val <- fits[[paste0(cohort_val, "_old")]][[graph_type]]
  #}

  # compare with largest connected component only
  #g2 <- get_largest_cc(g2)
  #if (!sentinel %in% graph::nodes(g2)) {
  #  g2 <- graph::addNode(sentinel, g2)
  #}

  # get adjacency matrices
  g_adj <- as(graph, "matrix")
  g_val_adj <- as(graph_val, "matrix")

  # ensure that we have the same nodes only in all graphs.
  # this might somewhat change results, but otherwise we
  # cant compute the MCC properly.
  use <- intersect(colnames(g_adj), colnames(g_val_adj))
  if (length(use) > 1) {
    g_adj <- g_adj[use, use]
    g_val_adj <- g_val_adj[use, use]

    # calculate performance using the DBgraph method compare()
    comp <- BDgraph::compare(g_adj, g_val_adj)
    f1 <- comp["F1-score", "estimate1"]
    mcc <- comp["MCC", "estimate1"]

    print(paste0("MCC: ", format(mcc, digits = 3)))
    print(paste0("F1: ", format(f1, digits = 3)))

    # the fraction of nodes retained in the overlap w.r.t. to the
    # total number of possible nodes
    common_nodes <- ncol(g_adj)
    mcc_frac <- common_nodes / ncol(data_val)
    row <- c(row, cross_cohort_f1=f1, cross_cohort_mcc=mcc,
             cross_cohort_mcc_frac=mcc_frac,
             common_nodes=common_nodes)
  } else {
    row <- c(row, cross_cohort_f1=NA, cross_cohort_mcc=NA,
             cross_cohort_mcc_frac=NA,
             common_nodes=NA)
  }

  # ------------------------------------------------------------------------------
  print("Checking mediation.")
  # ------------------------------------------------------------------------------
  if ("kora" %in% cohort) {
    med_val <- readRDS(fmediation_lolipop)
    med_fit <- readRDS(fmediation_kora)
  } else {
    med_val <- readRDS(fmediation_kora)
    med_fit <- readRDS(fmediation_lolipop)
  }
  row <-
    c(row,
      mediation.summary(med_val, sgenes, sgenes_selected, mediation_cutoff))

  # we also check the correspondence of the correlation values for all genes
  med_comparison <- compare_mediation_results(
    sentinel,
    med_val,
    med_fit,
    sgenes_selected,
    mediation_cutoff
  )

  row <- c(
    row, med_comparison
  )

  # ------------------------------------------------------------------------------
  print("Validating cis-eQTLs.")
  # ------------------------------------------------------------------------------

  # filter ceqtl to be only related to our sentinel SNP
  # TODO use proxy/high ld snps to increase ceqtl number?
  ceqtlsub <- ceqtl[ceqtl$Rs_ID %in% snp]
  if (nrow(ceqtlsub) < 1) {
    warning("Sentinel " %+% sentinel %+% " not found in cis-eQTL data")
    # report NA in stats file
    row <- c(row, cisEqtl=NA)
  } else {
    ceqtl_sgenes <- sgenes[sgenes %in% unlist(strsplit(ceqtlsub$Transcript_GeneSymbol, "\\|"))]
    ceqtl_sgenes_selected <-
      intersect(ceqtl_sgenes, sgenes_selected)

    # create matrix for fisher test
    cont <-
      matrix(
        c(
          length(ceqtl_sgenes),
          length(ceqtl_sgenes_selected),
          length(sgenes),
          length(sgenes_selected)
        ),
        nrow = 2,
        ncol = 2,
        byrow = T
      )
    rownames(cont) <- c("ceqtl", "no ceqtl")
    colnames(cont) <- c("not selected", "selected")
    row <- c(row,
             cisEqtl=fisher.test(cont)$p.value)
  }

  # ------------------------------------------------------------------------------
  print("CpG-gene and TF-CpG validation.")
  # ------------------------------------------------------------------------------
  # filter teqtl to be only related to our sentinel SNP
  # TODO use proxy/high ld snps to increase teqtl number?
  teqtlsub <-
    teqtl[teqtl$Rs_ID %in% snp]# & (teqtl$log10FDR) < (-2),,drop=F]
  if (nrow(teqtlsub) < 1) {
    warning("Sentinel " %+% sentinel %+% " not available in trans-eQTL data.")
    # report NA in stats file
    row <- c(row, transEqtl_tgenes=NA, transEqtl_tfs=NA)
  } else {
    if(ranges$seed == "meqtl") {
      row <- c(row,
               validate_trans_genes(teqtlsub, cgenes, tfs,
                                    cgenes_selected, tfs_selected))
    } else {
      row <- c(row,
               validate_trans_genes(teqtlsub, trans_entities, tfs,
                                    trans_entities_selected, tfs_selected))
    }
  }

  # we can also count how many of the identified cpg
  # genes are eQTMs in the bonder dataset
  # first convert bonder ensembl ids to symbols
  # (bonder results already filtered for sign. assoc.)
  bnd <- eqtms[, c("SNPName", "HGNCName")]
  bnd_symbol <- unique(bnd$HGNCName)
  row <-
    c(row, validate_cpggenes(bnd_symbol, cgenes, cgenes_selected))

  # ------------------------------------------------------------------------------
  print("Gene-Gene validation.")
  # ------------------------------------------------------------------------------
  # load geo whole blood expression data
  geo_sym <- all_genes[all_genes %in% geo$symbol]
  geosub <- geo[geo_sym]
  geosub <- t(geosub[, -1, with = F])
  colnames(geosub) <- geo_sym

  # list of all expr datasets
  expr.data <- list(geo = geosub, cohort = data_val[, all_genes, drop = F])
  val_g2g <- validate_gene2gene(expr.data, graph_maxcluster, all_genes)
  names(val_g2g) <- c("geo_gene_gene", "cohort_gene_gene")
  row <- c(row, val_g2g)

  # ----------------------------------------------------------------------------
  print("GO enrichment.")
  # ----------------------------------------------------------------------------
  # we want only genes which were selected by the model
  network_genes <- unique(c(cgenes_selected, tfs_selected,
                            sgenes_selected, spath_selected))
  if(ranges$seed != "meqtl") {
    network_genes <- unique(c(network_genes, trans_entities_selected))
  }

  go <- validate_geneenrichment(network_genes)
  if(!is.null(go)) {
    names(go) <- c("go_ids", "go_terms", "go_pvals", "go_qvals")
  } else {
    go <- rep(NA, 4)
    names(go) <- c("go_ids", "go_terms", "go_pvals", "go_qvals")
  }
  row <- c(row, go)

  return(row)
}, mc.cores=threads)

# write output file
valid <- do.call(rbind, valid)
if(is.null(dim(valid)) || ncol(valid) < 2) {
  stop("Error: wrong result dimensions.")
}

write.table(file=fout, valid, col.names=T,
            row.names=F, quote=F, sep="\t")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
