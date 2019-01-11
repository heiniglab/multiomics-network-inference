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
feqtl <- snakemake@input[[1]]
fkora_data <- snakemake@input$kora_data
flolipop_data <- snakemake@input$lolipop_data

# output
fout_plot <- snakemake@output$plot
dout_loci <- snakemake@output$loci_dir

# params

# minimum number of trans loci
hots_thres <- as.numeric(snakemake@wildcards$hots_thres)
reduce_region <- 10000

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
load(fkora_data)
available_snps <- colnames(geno)
#load(flolipop_data)
#available_snps <- intersect(available_snps, colnames(geno))

rm(geno,expr,meth,covars)
gc()

eqtl <- fread(paste0("zcat ", feqtl))

# filter eqtl for genes in our gene annotation
ga <- get.gene.annotation()
eqtl <- eqtl[eqtl$GeneSymbol %in% ga$SYMBOL,]
eqtl <- eqtl[eqtl$SNP %in% available_snps,]

# get all trans genes per SNP
trans_genes_by_snp <- tapply(eqtl$GeneSymbol, eqtl$SNP, function(x) {
  if(length(unique(x)) >= hots_thres) {
    return(unique(x))
  } else{
    return(NULL)
  }
})

# remove NULLs
trans_genes_by_snp <- trans_genes_by_snp[!unlist(lapply(trans_genes_by_snp,
                                                        is.null))]

# extract dataframe
hotspots <- cbind.data.frame(sentinel=names(trans_genes_by_snp),
                             ntrans=unlist(lapply(trans_genes_by_snp, length)),
                             stringsAsFactors=F)
hotspots <- cbind(hotspots,
                  eqtl[match(hotspots$sentinel, eqtl$SNP),
                       c("SNPPos", "SNPChr")])
hotspot_ranges <- with(hotspots, GRanges(SNPChr, IRanges(SNPPos, width=1)))
hotspot_ranges_reduced <- reduce(hotspot_ranges, min.gap=reduce_region,
                                 with.rev=T)

print("Total number of hotspots:")
print(nrow(hotspots))

print("Total number of reduced hotspots:")
print(length(hotspot_ranges_reduced))

# create new reduced hotspot list according to the reduced ranges
ntrans_reduced <- unlist(lapply(hotspot_ranges_reduced, function(r) {
  idxs <- unlist(r$revmap)
  spots <- hotspots[idxs,]
  ntrans <- length(unique(unlist(eqtl[eqtl$SNP %in% spots$sentinel,
                                      "GeneSymbol"])))
  ntrans
}))
ntrans_reduced <- cbind.data.frame(ntrans = ntrans_reduced)

# ------------------------------------------------------------------------------
print("Saving and plotting results.")
# ------------------------------------------------------------------------------
# create dummy files (more convenient for snakemake) for each sentinel
for(i in 1:nrow(hotspots)) {
  file.create(paste0(dout_loci, hotspots[i,"sentinel"], ".dmy"))
}

# plot a simple histogram for now
theme_set(theme_bw())

pdf(fout_plot)
ggplot(aes(x=ntrans), data=hotspots) + geom_histogram() +
  ggtitle(paste0("Overview on number of trans genes for ", nrow(hotspots), " hotspots.")) +
  xlab("number of trans associations")

ggplot(aes(x=ntrans), data=ntrans_reduced) + geom_histogram() +
  ggtitle(paste0("Overview on number of trans genes for ",
                 nrow(ntrans_reduced), " reduced hotspots (", reduce_region, "bp).")) +
  xlab("number of trans associations")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
