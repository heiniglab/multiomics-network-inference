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
library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(cowplot)

source("scripts/lib.R")

cols <- get_defaultcolors()
sfm <- scale_fill_manual(values=cols)
theme_set(theme_cowplot())

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
# input
finputs <- snakemake@input

# check seed group: eqtl or meqtl based
seeds <- ifelse(grepl("eqtlgen", finputs), "eqtl", "meqtl")

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
    c(sentinel, sg, tfs, sp, trans_assoc, cpg_genes, total, seed)
  } else {
    c(sentinel, sg, tfs, sp, trans_assoc, 0, total, seed)
  }

})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))

colnames(data) <- c("snp","snp_genes","TFs","shortest_path",
                      "trans_entities", "cpg_genes", "total", "seed")

# convert back to numeric..
data$snp_genes <- as.numeric(data$snp_genes)
data$trans_entities <- as.numeric(data$trans_entities)
data$TFs <- as.numeric(data$TFs)
data$shortest_path <- as.numeric(data$shortest_path)
data$cpg_genes <- as.numeric(data$cpg_genes)
data$total <- as.numeric(data$total)

# get fractions as well
data_fractions <- data
data_fractions$snp_genes <- data$snp_genes / data$total
data_fractions$trans_entities <- data$trans_entities / data$total
data_fractions$TFs <- data$TFs / data$total
data_fractions$shortest_path <- data$shortest_path / data$total
data_fractions$cpg_genes <- data$cpg_genes / data$total

# remove totals
data$total <- NULL
data_fractions$total <- NULL

# melt data frames for use in ggplot
melted <- melt(data)
melted_fractions <- melt(data_fractions)

# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------

gt <- ggtitle(paste0("Overview on number of entities gathered for\n",
                     length(finputs), " trans loci."))

fw <- facet_wrap( ~ seed, ncol=2)
vert_labels <- theme(axis.text.x = element_text(angle = 90, hjust = 1))

# violin plots containing points/lines showing distributions
gp <- ggplot(aes(y=value, x=variable, fill=variable), data=melted) +
	geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
	sfm + gt + fw + vert_labels
gp1 <- ggplot(aes(y=value, x=variable, fill=variable), data=melted_fractions) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  sfm + gt + fw + vert_labels
gp2 <- gp + geom_line(aes(group=snp)) + fw
# histograms of individual entity types
gp3 <- ggplot(aes(x=value, fill=variable), data=melted) +
  geom_histogram(stat="count") + facet_grid(seed ~ variable) + sfm

# barplot over all loci
gp4 <- ggplot(aes(y=value, x=snp, fill=variable), data=melted) +
  geom_bar(stat="identity", position="dodge") + sfm +
  vert_labels + fw

pdf(fout, width=10, height=8)
gp
gp1
gp2
gp3
gp4

# finally, plot the number of entities per locus
# TODO there must be a better way for this...
eqtl <- melted[melted$seed == "eqtl",]
meqtl <- melted[melted$seed == "meqtl",]
sum_per_locus_meqtl <- tapply(meqtl$value, meqtl$snp, sum)
sum_per_locus_meqtl <- cbind.data.frame(snp=names(sum_per_locus_meqtl),
                                       count=sum_per_locus_meqtl,
                                       stringsAsFactors=F)
sum_per_locus_eqtl <- tapply(eqtl$value, eqtl$snp, sum)
sum_per_locus_eqtl <- cbind.data.frame(snp=names(sum_per_locus_eqtl),
                                       count=sum_per_locus_eqtl,
                                       stringsAsFactors=F)

sum_per_locus <- rbind(sum_per_locus_meqtl, sum_per_locus_eqtl)
sum_per_locus$seed <- c(rep("meqtl", nrow(sum_per_locus_meqtl)),
                        rep("eqtl", nrow(sum_per_locus_eqtl)))
gp5 <- ggplot(aes(y=count, x="all loci", fill="all loci"),
              data = sum_per_locus) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("") + ylab("Number of entities") +
  ggtitle("Total number of entities for all available loci.") +
  scale_fill_manual(values=cols, guide=F) + facet_wrap( ~ seed, ncol=2)
gp5

dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
