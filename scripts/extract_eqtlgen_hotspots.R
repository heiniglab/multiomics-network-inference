#' -----------------------------------------------------------------------------
#' Extract eQTLgen hotspots from the eQTLgen trans eQTL file
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Fri Dec 21 21:40:48 2018
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))
library(data.table)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------
# input
feqtl <- snakemake@input$eqtl
fkora_data <- snakemake@input$kora_data
flolipop_data <- snakemake@input$lolipop_data
fgene_annot <- snakemake@input$gene_annot

# output
fout_plot <- snakemake@output$plot
dout_loci <- snakemake@output$loci_dir
fout_table <- snakemake@output$table

# params

# minimum number of trans loci
hots_thres <- as.numeric(snakemake@wildcards$hots_thres)

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
load(fkora_data)
available_snps <- colnames(geno)
load(flolipop_data)
available_snps <- intersect(available_snps, colnames(geno))

rm(geno,expr,meth,covars)
gc()

eqtl <- readRDS(feqtl)

# filter eqtl for genes in our gene annotation
ga <- load_gene_annotation(fgene_annot)
eqtl <- eqtl[eqtl$GeneSymbol %in% ga$SYMBOL,]
eqtl <- eqtl[eqtl$SNP %in% available_snps,]

# get all trans genes per SNP
trans_genes_by_snp <- tapply(eqtl$GeneSymbol, eqtl$SNP, function(x) {
  return(length(unique(x)))
})

ntrans <- unlist(trans_genes_by_snp, length)
eqtl$ntrans <- ntrans[match(eqtl$SNP, names(ntrans))]

# extract dataframe
hotspots <- eqtl[ntrans >= hots_thres]
hotspot_ranges <- unique(with(hotspots, GRanges(SNPChr, IRanges(SNPPos, width=1))))

print("Total number of hotspots:")
print(length(hotspot_ranges))

# ------------------------------------------------------------------------------
print("Saving and plotting results.")
# ------------------------------------------------------------------------------
fwrite(file=fout_table, hotspots, sep="\t", col.names=T)

# create dummy files (more convenient for snakemake) for each sentinel
for(i in 1:nrow(hotspots)) {
  file.create(paste0(dout_loci, hotspots[i,"SNP"], ".dmy"))
}

# plot a simple histogram for now
theme_set(theme_linedraw())
toplot <- hotspots[!duplicated(SNP)]

pdf(fout_plot)
ggplot(aes(x=ntrans), data=hotspots) + geom_histogram() +
  ggtitle(paste0("Overview on number of trans genes for ", nrow(hotspots), " hotspots.")) +
  xlab("number of trans associations")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
