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
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(Homo.sapiens))
suppressPackageStartupMessages(library(illuminaHumanv3.db))
suppressPackageStartupMessages(library(ggplot2))
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
graph_type <- snakemake@wildcards$graph_type
cohort <- snakemake@wildcards$cohort

# output
# main outfile for validation results
fout <- snakemake@output[[1]]

# additional outputs, mostly plots
fmediation_summary_plot <- snakemake@output$mediation_summary_plot
fgraphdot <- snakemake@output[["graphdot"]]
fgraphpdf <- snakemake@output[["graphpdf"]]

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------

print("Loading gene expression data.")

# load GEO/ARCHS4  data
geo <- fread(fgeo,
             header = T, sep = "\t")
colnames(geo)[1] <- "symbol"
setkey(geo, symbol)

#print("Loading KORA eQTL.")
#kora_ceqtl <- fread(fciseqtl_kora, sep = "\t",
#                    header = T)
#kora_teqtl <- fread(ftranseqtl_kora, sep = "\t",
#                    header = T)

print("Loading Joehanes eQTL.")
# load cis-eQTL data
ceqtl <- fread(paste0("zcat ", fciseqtl_joehanes),
               sep = ",")
# get only cis eqtls defined in the paper
ceqtl <- ceqtl[ceqtl$Is_Cis == 1]
print(paste0("Loaded ", nrow(ceqtl), " cis-eQTL"))

# load the trans-eQTL dataset
teqtl <- fread(paste0("zcat ", ftranseqtl_joehanes),
               sep = ",")

# sanity check
print(paste0("Loaded ", nrow(teqtl), " trans-eQTL"))

# load the bonder cis-eQTMs for cpg-gene validation
# TODO we could use those for prior definitions, then use the ROADMAP
# chromhmm states for validation...
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

ranges <- readRDS(franges)

# ------------------------------------------------------------------------------
print("Loading GGM fits.")
# ------------------------------------------------------------------------------
kfit <- readRDS(fkora_fit)
lfit <- readRDS(flolipop_fit)
fits <- list(kora = kfit, lolipop = lfit)

print(paste0("Validating ", cohort, " fit for '", graph_type , "' graph fit."))

row <- c(sentinel = sentinel,
         cohort = cohort,
         graph_type = graph_type)

# ------------------------------------------------------------------------------
print("Preparing fit and validation data.")
# ------------------------------------------------------------------------------
g <- fits[[cohort]][[graph_type]]
ds <- graph::degree(g)
keep <- names(ds[ds>0])
g <- subGraph(keep, g)

print("Plotting graph.")
pdf(fgraphpdf)
temp <- plot.ggm(g, sentinel, T, dot.out = fgraphdot)
dev.off()

# we get the according dataset on which to validate
if ("lolipop" %in% cohort) {
  # ggm calculated on lolipop, validate on kora
  data_val <- kdata
  data_fit <- ldata
} else {
  # assume kora cohort, validate on lolipop
  data_val <- ldata
  data_fit <- kdata
}
# dnodes -> full set of possible nodes
dnodes <- colnames(data_val)

# ------------------------------------------------------------------------------
print("Getting basic stats (number nodes, number edges, graph density)")
# ------------------------------------------------------------------------------
nn <- numNodes(g)
ne <- numEdges(g)
gd <- (ne * 2) / (nn * (nn - 1))
row <- c(row, number_nodes=nn, number_edges=ne, graph_density=gd)

# ------------------------------------------------------------------------------
print("Getting cluster information")
# ------------------------------------------------------------------------------
ig = igraph::graph_from_graphnel(g)
cl = clusters(ig)
ncluster <- cl$no
scluster <- paste(cl$csize, collapse = ",")

# remember snp membership
snp_cluster <- NA
if (sentinel %in% names(cl$membership)) {
  snp_cluster <- cl$membership[sentinel]
  snp_cluster_size <- cl$csize[snp_cluster]
  snp_cluster <- paste0(snp_cluster, "-", snp_cluster_size)
}
row <- c(row,
         cluster=ncluster, cluster_sizes=scluster, snp_cluster=snp_cluster)

# ------------------------------------------------------------------------------
print("Getting largest CC for validation.")
# ------------------------------------------------------------------------------
keep <- cl$membership == which.max(cl$csize)
keep <- names(cl$membership[keep])
g <- graph::subGraph(keep, g)

# the nodes retained in the fitted graph model in the largest CC
gnodes <- graph::nodes(g)
if (!sentinel %in% gnodes) {
  g <- graph::addNode(sentinel, g)
}

# ------------------------------------------------------------------------------
print("Calculating graph score.")
# ------------------------------------------------------------------------------
# we use the (full) igraph object for this, will be filtered for the sentinel
# cluster
score <- get_graph_score(ig, sentinel, ranges)
row <- c(row,
         graph_score = score)

# ------------------------------------------------------------------------------
print("Defining entity sets (selected / not selected)")
# ------------------------------------------------------------------------------

# for now, filter only for those nodes for which data was available
# should change once we update the data matrix...
gnodes <- gnodes[gnodes %in% dnodes]

# get names of all entities, total and selected by ggm
snp <- sentinel
data_val[, snp] <- as.integer(as.character(data_val[, snp]))
data_fit[, snp] <- as.integer(as.character(data_fit[, snp]))
if(ranges$seed == "meqtl") {
  trans_entities <- intersect(dnodes, names(ranges$cpgs))
} else {
  trans_entities <- intersect(dnodes,  ranges$trans_genes$SYMBOL)
}
trans_entities_selected <- trans_entities[trans_entities %in% gnodes]

all_genes <- dnodes[!grepl("^rs|^cg", dnodes)]
sgenes <- intersect(dnodes, ranges$snp_genes$SYMBOL)
if (!is.na(snp_cluster)) {
  sgenes_selected <- sgenes[sgenes %in% unlist(adj(g, snp))]
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
    cgenes_selected <- cgenes[cgenes %in% unlist(adj(g, trans_entities_selected))]
  } else {
    cgenes_selected <- c()
  }
} else {
  cgenes_selected <- c()
}

if(length(trans_entities_selected) > 0) {
  tfs_selected <- tfs[tfs %in% unlist(adj(g, trans_entities_selected))]
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
  spath_selected=length(spath_selected)
)

# ------------------------------------------------------------------------------
print("Using MCC to check how well graph replicated across cohorts.")
# ------------------------------------------------------------------------------
if ("lolipop" %in% cohort) {
  cohort2 <- "kora"
} else {
  cohort2 <- "lolipop"
}

# get graph fit on other cohort
g2 <- fits[[cohort2]][[graph_type]]

# compare with largest connected component only
g2 <- get_largest_cc(g2)
if (!sentinel %in% graph::nodes(g2)) {
  g2 <- graph::addNode(sentinel, g2)
}

# get adjacency matrices
g_adj <- as(g, "matrix")
g2_adj <- as(g2, "matrix")

# ensure that we have the same nodes only in all graphs.
# this might somewhat change results, but otherwise we
# cant compute the MCC properly.
use <- intersect(colnames(g_adj), colnames(g2_adj))
if (length(use) > 1) {
  g_adj <- g_adj[use, use]
  g2_adj <- g2_adj[use, use]

  # calculate performance using the DBgraph method compare()
  mcc <- BDgraph::compare(g_adj, g2_adj)["MCC", "estimate"]
  print(paste0("MCC: ", format(mcc, digits = 3)))

  # the fraction of nodes retained in the overlap w.r.t. to the
  # total number of possible nodes
  mcc_frac <- ncol(g_adj) / ncol(data_val)
  row <- c(row, cross_cohort_mcc=mcc, cross_cohort_mcc_frac=mcc_frac)
} else {
  row <- c(row, cross_cohort_mcc=NA, cross_cohort_mcc_frac=NA)
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
row <- c(
  row,
  compare_mediation_results(
    sentinel,
    med_val,
    med_fit,
    sgenes_selected,
    mediation_cutoff,
    fmediation_summary_plot
  )
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
  row <- c(row, transEqtl_cgenes=NA, transEqtl_tfs=NA)
} else {
  row <- c(row,
           validate_trans_genes(teqtlsub, cgenes, tfs,
                                cgenes_selected, tfs_selected))
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
val_g2g <- validate_gene2gene(expr.data, g, all_genes)
names(val_g2g) <- c("geo_gene_gene", "cohort_gene_gene")
row <- c(row, val_g2g)

# ------------------------------------------------------------------------------
print("GO enrichment.")
# ------------------------------------------------------------------------------
go <- validate_geneenrichment(gnodes)
if(!is.null(go)) {
  names(go) <- c("go_ids", "go_terms", "go_pvals", "go_qvals")
} else {
  go <- rep(NA, 4)
  names(go) <- c("go_ids", "go_terms", "go_pvals", "go_qvals")
}
row <- c(row, go)

# write output file
write.table(file=fout, t(row), col.names=T,
            row.names=F, quote=F, sep="\t")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
