# ------------------------------------------------------------------------------
#'
#' Script to collect all associated entities (expr and geno)
#' for a specific sentinel snp of the GTEx trans- eQTL results
#'
#' @author Johann Hawe
#'
#' @date 2017/03/04
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load packages, source scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
library(data.table)
library(graph)
source("scripts/lib.R")
source("scripts/collect_ranges_methods.R")

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
# inputs
feqtl <- snakemake@input$eqtl
fppi <- snakemake@input$ppi
ftfbs_annot <- snakemake@input$tfbs_annot
fgene_annot <- snakemake@input$gene_annot

# output file
fout_ranges <- snakemake@output$ranges
fout_plot <- snakemake@output$plot

# params
sentinel <- snakemake@wildcards$sentinel
tissue <- snakemake@wildcards$tissue

# ------------------------------------------------------------------------------
print("Load and preprocess data.")
# ------------------------------------------------------------------------------

# PPIs
ppi_db <- readRDS(fppi)
ppi_genes <- nodes(ppi_db)

# gene annotation
gene_annot <- load_gene_annotation(fgene_annot)
gene_annot$ids <- probes.from.symbols(gene_annot$SYMBOL,
                                      as.list=T)

# load trans-eQTL
eqtl = fread(feqtl)
eqtl <- eqtl[eqtl$SNP == sentinel,]
if(nrow(eqtl)<5) stop("Will only process hotspots with at least 5 associations.")

# load TFBS annotation
tfbs <- readRDS(ftfbs_annot)

# ------------------------------------------------------------------------------
print("Collecting SNP and gene ranges.")
# ------------------------------------------------------------------------------

# get sentinel region + extended region in which to look for genes
spos <- eqtl[1, c("SNPChr", "SNPPos")]
if(nrow(spos)<1) {
  stop("Couldn't get SNP pos.")
}

sentinel_range <- with(spos,
                       GRanges(paste0("chr", SNPChr),
                               IRanges(SNPPos, width=1)))
sentinel_extrange <- resize(sentinel_range, 1e6, fix="center")
names(sentinel_range) <- names(sentinel_extrange) <- sentinel

# get the regions of the associated trans genes
# we dont have the end position/length of the respective genes in the eqtl
# table, so we use our own annotation
trans_genes <- eqtl$GeneSymbol
trans_genes <- gene_annot[gene_annot$SYMBOL %in% trans_genes]

# could be some are missing. in that case we very likely cannot do anything
# aboutit, since we also won't have any probe ids. if there are too few left,
# we have to abort and report an error
if(length(trans_genes) < 5) {
  stop("Too many trans genes missing in our annotation. Will not collect
       ranges.")
}

# ------------------------------------------------------------------------------
print("Retrieving SNP genes.")
# ------------------------------------------------------------------------------

# get the relevant snp genes by overlapping with our sentinel region
snp_genes <- subsetByOverlaps(promoters(gene_annot),
                               sentinel_extrange,
                               ignore.strand=T)

# ------------------------------------------------------------------------------
print("Collecting TFs and shortest path genes.")
# ------------------------------------------------------------------------------

# load all TFBS we have available in our data and connect with trans-genes
tfbs <- tfbs[rownames(tfbs) %in% names(trans_genes),,drop=F]
tfs <- NULL
sp <- NULL

tfs_by_transGene <- get_tfs_by_transGene(tfbs, trans_genes, gene_annot)
if(length(tfs_by_transGene) > 0) {
  tfs <- unique(unlist(GenomicRangesList(tfs_by_transGene), use.names=F))
}

# find the shortest path genes between the SNP genes and the annotated TFs
if(length(tfs)<1){
  warning("No TFs, skipping shortest paths calculation.")
} else {
  sp <- collect_shortest_path_genes(tfs$SYMBOL, trans_genes$SYMBOL,
                                    tfs_by_transGene, ppi_genes, 
                                    snp_genes$SYMBOL, ppi_db, gene_annot)
}

# ------------------------------------------------------------------------------
print("Finalizing and saving results.")
# ------------------------------------------------------------------------------
result <-  list(sentinel = sentinel_range,
                snp_genes=snp_genes,
                trans_genes=trans_genes)

if(!is.null(sp)){
  result$spath <- sp
}
if(!is.null(tfs)){
  result$tfs <- tfs
  result$tfs_by_transGene <- tfs_by_transGene
}

# set seed name
result$seed <- "eqtl"
result$tissue <- tissue

saveRDS(file=fout_ranges, result)

# plot
plot_ranges(result, fout_plot)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
