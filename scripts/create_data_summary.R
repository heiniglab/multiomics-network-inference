# ------------------------------------------------------------------------------
#' Script to create a summary plot over the data of all loci.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Loading libraries and sourcing scripts.")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
library(ggplot2)
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
  # load data
  data <- readRDS(fin)

  # for now just report the number of identified entities
  ncol(data)
})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))
colnames(data) <- c("entities")
data$seed <- ifelse(grepl("eqtlgen", finputs), "eqtlgen", "meqtl")

# convert back to numeric..
data$entities <- as.numeric(data$entities)

# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------
pdf(fout)
ggplot(data, aes(x=entities)) + geom_histogram(stat="count") +
  xlab("number of variables") + ggtitle("Histogram of the number of
                                                     variables over all loci
                                                     for meQTL and eQTLgen.") +
  facet_wrap(~seed, ncol=2)
dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
