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

# check seed group: twas or eqtl based
group <- if(grepl("twas", finputs), "twas", "eqtl")

# output
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
  tfs <- length(ranges$tfs)
  sp <- length(ranges$spath)

  seed <- ranges$seed

  if(seed == "meqtl"){
    trans_assoc <- length(ranges$cpgs)
    cpg_genes <- length(ranges$cpg_genes)
  } else {
    trans_assoc <- length(ranges$trans_genes)
    cpg_genes <- NA
  }

  # get the grand total of entities
  total <- sg + tfs + sp + trans_assoc + ifelse(is.na(cpg_genes), 0, cpg_genes)

  if(seed == "meqtl") {
    c(sentinel, sg, tfs, sp, trans_assoc, cpg_genes, total)
  } else {
    c(sentinel, sg, tfs, sp, trans_assoc, total)
  }

})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))

if(seed == "meqtl") {
  colnames(data) <- c("snp","snp_genes","TFs","shortest_path",
                      "trans_entities", "cpg_genes", "total")
} else {
  colnames(data) <- c("snp","snp_genes","TFs","shortest_path", "trans_entities",
                      "total")
}

# convert back to numeric..
data$snp_genes <- as.numeric(data$snp_genes)
data$trans_entities <- as.numeric(data$trans_entities)
data$TFs <- as.numeric(data$TFs)
data$shortest_path <- as.numeric(data$shortest_path)
if(seed == "meqtl") data$cpg_genes <- as.numeric(data$cpg_genes)
data$total <- as.numeric(data$total)

# get fractions as well
data_fractions <- data
data_fractions$snp_genes <- data$snp_genes / data$total
data_fractions$trans_entities <- data$trans_entities / data$total
data_fractions$TFs <- data$TFs / data$total
data_fractions$shortest_path <- data$shortest_path / data$total
if(seed == "meqtl") data_fractions$cpg_genes <- data$cpg_genes / data$total

# remove totals
data$total <- NULL
data_fractions$total <- NULL

# melt the data frame for use in ggplot
melted <- melt(data)
colnames(melted) <- c("locus", "entity", "count")
melted_fractions <- melt(data_fractions)
colnames(melted_fractions) <- c("locus", "entity", "fraction")


# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------

gt <- ggtitle(paste0("Overview on number of entities gathered for\n",
                     length(finputs), " trans loci."))
# violin plots containing points/lines showing distributions
gp <- ggplot(aes(y=count, x=entity, fill=entity), data=melted) +
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	sfm + gt
gp1 <- ggplot(aes(y=fraction, x=entity, fill=entity), data=melted_fractions) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  sfm + gt
gp2 <- gp + geom_line(aes(group=locus))
# histograms of individual entity types
gp3 <- ggplot(aes(x=count, fill=entity), data=melted) +
  geom_histogram(stat="count") + facet_wrap(~ entity, ncol=2) + sfm

# barplot over all loci
gp4 <- ggplot(aes(y=count, x=locus, fill=entity), data=melted) +
  geom_bar(stat="identity", position="dodge") + sfm +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf(fout, width=10, height=8)
gp
gp1
gp2
gp3
gp4

# finally, plot the number of entities per locus
sum_per_locus <- tapply(melted$count, melted$locus, sum)
sum_per_locus <- cbind.data.frame(locus=names(sum_per_locus),
                       count=sum_per_locus, stringsAsFactors=F)
gp5 <- ggplot(aes(y=count, x="all loci", fill="all loci"),data = sum_per_locus) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("") + ylab("Number of entities") +
  ggtitle("Total number of entities for all available loci.") +
  scale_fill_manual(values=cols, guide=F)
gp5

dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
