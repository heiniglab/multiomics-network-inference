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
fmeqtl <- snakemake@input[[1]]

# output
fout_plot <- snakemake@output$plot
dout_loci <- snakemake@output$loci_dir

# params
threads <- snakemake@threads

# minimum number of trans loci
hots_thres <- as.numeric(snakemake@wildcards$hots_thres)

# ------------------------------------------------------------------------------
print("Loading and processing data.")
# ------------------------------------------------------------------------------
meqtl <- fread(fmeqtl)

# check the number of trans sentinel cpg for each sentinel
trans_cpgs_by_snp <- tapply(meqtl$sentinel.cpg, meqtl$sentinel.snp, function(x){
  if(length(unique(x)) >= hots_thres) {
    return(unique(x))
  } else {
    return(NULL)
  }
})

# remove NULLs
trans_cpgs_by_snp <- trans_cpgs_by_snp[!unlist(lapply(trans_cpgs_by_snp,
                                                      is.null))]

# extract dataframe
hotspots <- cbind.data.frame(sentinel=names(trans_cpgs_by_snp),
                             ntrans=unlist(lapply(trans_cpgs_by_snp, length)),
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
hist(hotspots$ntrans, breaks=50, xlab="number of trans genes",
     main = "Hotspot overview")
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
