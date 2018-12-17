# ------------------------------------------------------------------------------
#' Script to create a summary plot over all generated loci.
#' Extracts all individual entities from the locus-ranges objects
#' and shows their numbers in an overview.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Loading libraries and sourcing scripts.")
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape)
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")
cols <- set_defaultcolors()
sfm <- scale_fill_manual(values=cols)
theme_set(theme_bw())

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
finputs <- snakemake@input
fout <- snakemake@output[[1]]

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
  tg <- length(ranges$trans_genes)
  tfs <- length(ranges$tfs)
  sp <- length(ranges$spath)

  c(sentinel, sg, tg, tfs, sp)
})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))
colnames(data) <- c("snp","snp_genes","trans_genes","TFs","shortest_path")

# convert back to numeric..
data$snp_genes <- as.numeric(data$snp_genes)
data$trans_genes <- as.numeric(data$trans_genes)
data$TFs <- as.numeric(data$TFs)
data$shortest_path <- as.numeric(data$shortest_path)

# melt the data frame for use in ggplot
melted <- melt(data)
colnames(melted) <- c("locus", "entity", "count")

# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------

gt <- ggtitle(paste0("Overview over all entities gathered for\n",
                     length(finputs), " trans-eQTL loci."))
# violin plots containing points/lines showing distributions
gp <- ggplot(aes(y=count, x=entity, fill=entity), data=melted) +
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	sfm + gt
gp1 <- gp + geom_jitter(aes(color=entity),height = 0, width = 0.1, size=0.1) +
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
dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
