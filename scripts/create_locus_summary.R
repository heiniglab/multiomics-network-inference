# ------------------------------------------------------------------------------
# Script to create a summary plot over all generated loci.
# Extracts all individual entities from the locus-ranges objects
# and shows their numbers in an overview.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Loading libraries and sourcing scripts.")
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape)
library(GenomicRanges)
source("scripts/lib.R")

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
  sentinel <- names(ranges$sentinel_range)
  # get number of entities
  cpgs <- length(ranges$cpgs)
  sg <- length(ranges$snp_genes)
  cg <- length(ranges$cpg_genes)
  tfs <- length(ranges$tfs)
  sp <- length(ranges$spath)

  c(sentinel, cpgs, sg, cg, tfs, sp)
})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))
colnames(data) <- c("snp", "cpgs","snp_genes","cpg_genes","TFs","shortest_path")

# convert back to numeric..
data$cpgs <- as.numeric(data$cpgs)
data$snp_genes <- as.numeric(data$snp_genes)
data$cpg_genes <- as.numeric(data$cpg_genes)
data$TFs <- as.numeric(data$TFs)
data$shortest_path <- as.numeric(data$shortest_path)

# melt the data frame for use in ggplot
melted <- melt(data)
colnames(melted) <- c("locus", "entity", "count")
# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------
gp <- ggplot(aes(y=count, x=entity, fill=entity), data=melted) + 
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
	geom_jitter(height = 0, width = 0.1) + 
	ggtitle(paste0("Overview over all entities gathered for\n", length(finputs), " network loci."))
ggsave(file=fout, gp, width=10, height=8)

