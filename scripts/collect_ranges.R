#'
#' Script to collect all associated ids (expr, meth and geno)
#' for a specific sentinel snp
#'
#' @author Johann Hawe
#'
#' @date 2017/03/04
#'

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
# load packages, source scripts
# ------------------------------------------------------------------------------
library(GenomicRanges)
library(GenomicFeatures)
library(FDb.InfiniumMethylation.hg19)
library(data.table)
library(illuminaHumanv3.db)
library(rtracklayer)
library(graph)
library(RBGL) # for shortest paths
library(Matrix)
source("scripts/lib.R")
source("scripts/collect_ranges_methods.R")

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------
fcosmo <- snakemake@input$tcosmo
fmeqtl <- snakemake@input$meqtl
fppi_db <- snakemake@input$ppi_db
fprio_tab <- snakemake@input$priorization
fgene_annot <- snakemake@input$gene_annot

# TODO: create this file from scratch!
fcpgcontext <- snakemake@input$cpgcontext

ofile <- snakemake@output[[1]]
sentinel <- snakemake@wildcards$sentinel

# ------------------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------------------

print("Loading data.")

gene_annot <- load_gene_annotation(fgene_annot)
gene_annot$ids <- probes.from.symbols(gene_annot$SYMBOL,
                                           as.list=T)
ppi_db <- readRDS(fppi_db)

# load trans-meQTL table
trans_meQTL = read.csv(fmeqtl,
                       sep="\t", stringsAsFactors=F)

# load trans-cosmo information
cosmo <- readRDS(fcosmo)


# load priorization table
prio <- read.table(fprio_tab, sep="\t", header=T, stringsAsFactors = F)

# get trans-genes which should be used for shortest path extraction
prio <- prio[prio$sentinel == sentinel,,drop=F]
if(nrow(prio) > 0) {
  best_trans <- unique(prio$trans.gene)
} else {
  best_trans <- NULL
}

# ------------------------------------------------------------------------------
# Collect and save ranges
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Collecting SNP and CpG ranges.")
# ------------------------------------------------------------------------------

pairs = which(trans_meQTL[,"sentinel.snp"] == sentinel)

# get sentinel idx
sentinel_idx <- which(cosmo$snp==sentinel)[1]

# the large interval range for the sentinel
chr <- paste0("chr", trans_meQTL[pairs,"chr.snp"][1])
start <- trans_meQTL[pairs,"interval.start.snp"][1]
end <- trans_meQTL[pairs,"interval.end.snp"][1]

sentinel_extrange <- GRanges(chr, IRanges(start,end))
sentinel_range <- with(cosmo[sentinel_idx,],
                       GRanges(paste0("chr", snp.chr),
                               IRanges(snp.pos, width=1)))
names(sentinel_range) <- names(sentinel_extrange) <- sentinel

# get cosmo subset
idxs <- get.trans.cpgs(sentinel, trans_meQTL, cosmo)

# get related genes, i.e. genes near meQTL loci (snp+cpg)

# extend cpgs
cosmosub <- cosmo[idxs,]

croi <- with(cosmosub, GRanges(paste0("chr", cpg.chr),
                               IRanges(cpg.pos,width=2)))
names(croi) <- as.character(cosmosub[,"cpg"])
croi <- unique(croi)

# extended sentinel region
sroi <- sentinel_extrange
names(sroi) <- sentinel

# ------------------------------------------------------------------------------
print("Retrieving SNP and CpG genes.")
# ------------------------------------------------------------------------------

# get the relevant snp genes by overlapping with our sentinel region
genes_sroi <- subsetByOverlaps(gene_annot, sroi, ignore.strand=T)

# get genes near our cpg regions
genes_by_cpg <- get.nearby.ranges(croi, promoters(gene_annot))
names(genes_by_cpg) <- names(croi)
# get original ranges (not promoters)
genes_by_cpg <- lapply(names(genes_by_cpg), function(cg) {
  gs <- genes_by_cpg[[cg]]
  gene_annot[gs$hit_idx]
})
names(genes_by_cpg) <- names(croi)

# get single list of all cpg genes
genes_croi <- unique(unlist(GRangesList(unlist(genes_by_cpg))))

# ------------------------------------------------------------------------------
print("Collecting TFs and shortest path genes.")
# ------------------------------------------------------------------------------

tfs <- NULL
sp <- NULL

# get cpg ids and SNP gene symbols
cpgs <- names(croi)
snp_genes <- unique(genes_sroi$SYMBOL)

# modify ppi_db to contain our CpGs

# load the cpg-tf context
tfbs_ann <- get_tfbs_context(names(croi), fcpgcontext)
cpgs_with_tfbs <- cpgs[cpgs %in% rownames(tfbs_ann[rowSums(tfbs_ann)>0,])]
snp_genes_in_string <- snp_genes[snp_genes %in% nodes(ppi_db)]

# get locus graph
locus_graph <- add.to.graphs(list(ppi_db), sentinel, snp_genes,
                             cpgs_with_tfbs, tfbs_ann)[[1]]

# get tfs connected to cpgs
tfs = unique(unlist(adj(locus_graph, cpgs_with_tfbs)))
print(paste0("Annotated TFs: ", paste(tfs, collapse=", ")))

if(length(tfs)<1){
  warning("No TFs, skipping shortest paths calculation.")
} else {
  # the nodes we want to keep
  nodeset <- c(nodes(ppi_db), setdiff(tfs, "KAP1"),
               snp_genes_in_string, cpgs_with_tfbs)
  locus_graph <- subGraph(intersect(nodes(locus_graph), nodeset), locus_graph)

  syms_sp <- get_shortest_paths(cis = cpgs_with_tfbs,
                                       trans=unique(c(snp_genes_in_string,
                                                      tfs)),
                                       snp_genes=snp_genes_in_string,
                                       locus_graph,
                                       best_trans)

  if(length(syms_sp) < 1){
    warning("No shortest path genes.")
  } else {
    sp <- gene_annot[gene_annot$SYMBOL %in% syms_sp]
  }
  tfs <- gene_annot[gene_annot$SYMBOL %in% tfs]
}

# ------------------------------------------------------------------------------
print("Finalizing and saving results.")
# ------------------------------------------------------------------------------
result <-  list(cpgs=croi,sentinel=sentinel_range,
                sentinel_ext_range=sentinel_extrange,
                snp_genes=genes_sroi, cpg_genes=genes_croi,
                cpg_genes_by_cpg=genes_by_cpg)

if(!is.null(sp)){
  result$spath <- sp
}
if(!is.null(tfs)){
  result$tfs <- tfs
}
# set seed name
result$seed <- "meqtl"

saveRDS(file=ofile, result)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
