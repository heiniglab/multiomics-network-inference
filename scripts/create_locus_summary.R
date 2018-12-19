# ------------------------------------------------------------------------------
#' Script to create a summary plot over all generated hotspot networks for a
#' specific seed.
#' Extracts all individual entities from the locus-ranges objects
#' and shows their numbers in an overview.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Loading libraries and sourcing scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
library(ggplot2)
library(reshape)
source("scripts/lib.R")

cols <- set_defaultcolors()
sfm <- scale_fill_manual(values=cols)
theme_set(theme_bw())

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
# input
finputs <- snakemake@input
# output
fout <- snakemake@output[[1]]

# read SEED type from first input file
seed <- readRDS(finputs[[1]])$seed

# ------------------------------------------------------------------------------
print("Loading data and creating data-frame for plotting.")
# ------------------------------------------------------------------------------
data <- lapply(finputs, function(fin) {
  # load ranges
  ranges <- readRDS(fin)

  # get sentinel name
  sentinel <- names(ranges$sentinel)
  # get number of entities
  sg <- length(ranges$snp_genes)
  tfs <- length(ranges$tfs)
  sp <- length(ranges$spath)

  if(seed == "meqtl"){
    trans_assoc <- length(ranges$cpgs)
    cpg_genes <- length(ranges$cpg_genes)
  } else {
    trans_assoc <- length(ranges$trans_genes)
    cpg_genes <- NULL
  }

  c(sentinel, sg, tfs, sp, trans_assoc, cpg_genes)
})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))

if(seed == "meqtl") {
  colnames(data) <- c("snp","snp_genes","TFs","shortest_path",
                      "trans_entities", "cpg_genes")
} else {
  colnames(data) <- c("snp","snp_genes","TFs","shortest_path", "trans_entities")
}

# convert back to numeric..
data$snp_genes <- as.numeric(data$snp_genes)
data$trans_entities <- as.numeric(data$trans_entities)
data$TFs <- as.numeric(data$TFs)
data$shortest_path <- as.numeric(data$shortest_path)
if(seed == "meqtl") data$cpg_genes <- as.numeric(data$cpg_genes)

# melt the data frame for use in ggplot
melted <- melt(data)
colnames(melted) <- c("locus", "entity", "count")

# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------

gt <- ggtitle(paste0("Overview over all entities gathered for\n",
                     length(finputs), " trans loci."))
# violin plots containing points/lines showing distributions
gp <- ggplot(aes(y=count, x=entity, fill=entity), data=melted) +
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	sfm + gt
gp1 <- gp + geom_jitter(aes(color=entity), height = 0, width = 0.1, size=0.1) +
  scale_color_manual(values=cols)
gp2 <- gp + geom_line(aes(group=locus))
# histograms of individual entity types
gp3 <- ggplot(aes(x=count, fill=entity), data=melted) +
  geom_histogram(stat="count") + facet_wrap(~ entity, ncol=2) + sfm
pdf(fout, width=10, height=8)
gp
gp1
gp2
gp3

# finally, plot the number of entities per locus
sum_per_locus <- tapply(melted$count, melted$locus, sum)
sum_per_locus <- cbind.data.frame(locus=names(sum_per_locus),
                       count=sum_per_locus, stringsAsFactors=F)
gp4 <- ggplot(aes(y=count, x="all loci", fill="all loci"),data = sum_per_locus) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("") + ylab("Number of entities") +
  ggtitle("Total number of entities for all available loci.") +
  scale_fill_manual(values=cols, guide=F)
gp4

dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
