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
fmeqtl <- snakemake@input[[1]]
fkora_data <- snakemake@input$kora_data
flolipop_data <- snakemake@input$lolipop_data

# output
fout_plot <- snakemake@output$plot
fout_table <- snakemake@output$table
dout_loci <- snakemake@output$loci_dir

# params
threads <- snakemake@threads

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

meqtl <- fread(fmeqtl)

# check the number of trans sentinel cpg for each sentinel
trans_cpgs_by_snp <- tapply(meqtl$sentinel.cpg, meqtl$sentinel.snp, function(x){
  return(length(unique(x)))
})

ntrans <- unlist(trans_cpgs_by_snp, length)
meqtl$ntrans <- ntrans[match(meqtl$sentinel.snp, names(ntrans))]

# extract dataframe
hotspots <- meqtl[ntrans >= hots_thres]
hotspots <- hotspots[hotspots$sentinel.snp %in% available_snps,,drop=F]

print("Total number of hotspots:")
print(length(unique(hotspots)))

# ------------------------------------------------------------------------------
print("Saving and plotting results.")
# ------------------------------------------------------------------------------
fwrite(file=fout_table, hotspots, sep="\t", col.names=T)

# create dummy files (more convenient for snakemake) for each sentinel
sents <- unique(hotspots$sentinel.snp)
for(i in 1:length(sents)) {
  file.create(paste0(dout_loci, sents[i], ".dmy"))
}


# plot a simple histogram for now
theme_set(theme_linedraw())
toplot <- hotspots[!duplicated(sentinel.snp)]

pdf(fout_plot)
ggplot(aes(x=ntrans), data=toplot) + geom_histogram() +
  ggtitle(paste0("Overview on number of trans CpGs for ", length(sents), " hotspots.")) +
  xlab("number of trans associations")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
