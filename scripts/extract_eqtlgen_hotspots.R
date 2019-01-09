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
library(data.table)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params")
# ------------------------------------------------------------------------------
# input
feqtl <- snakemake@input[[1]]
fkora_data <- snakemake@input$kora_data

# output
fout_plot <- snakemake@output$plot
dout_loci <- snakemake@output$loci_dir

# params

# minimum number of trans loci
hots_thres <- as.numeric(snakemake@wildcards$hots_thres)

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
load(fkora_data)
available_snps <- colnames(geno)
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
    return(x)
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
print("Total number of hotspots:")
print(nrow(hotspots))

# ------------------------------------------------------------------------------
print("Saving and plotting results.")
# ------------------------------------------------------------------------------
# create dummy files (more convenient for snakemake) for each sentinel
for(i in 1:nrow(hotspots)) {
  file.create(paste0(dout_loci, hotspots[i,"sentinel"], ".dmy"))
}

# plot a simple histogram for now
pdf(fout_plot)
hist(hotspots$ntrans, breaks=200, xlab="number of trans genes",
     main = "Hotspot overview")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
