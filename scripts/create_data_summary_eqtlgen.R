# ------------------------------------------------------------------------------
#' Script to create a summary plot over the data of all loci.
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
  # load data
  data <- readRDS(fin)

  # for now just report the number of identified entities
  ncol(data)
})

data <- do.call(rbind.data.frame, args=c(data, stringsAsFactors=F))
colnames(data) <- c("entities")

# convert back to numeric..
data$entities <- as.numeric(data$entities)

# ------------------------------------------------------------------------------
print("Plotting and saving results.")
# ------------------------------------------------------------------------------
pdf(fout)
hist(data$entities, breaks=100)
dev.off()

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
