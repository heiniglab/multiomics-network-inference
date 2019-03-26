# ------------------------------------------------------------------------------
#'
#' Script to validate the GGMs calculated on simulated
#' data
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------

sink(file=snakemake@log[[1]])

# ------------------------------------------------------------------------------
# load needed libraries and source scripts
# ------------------------------------------------------------------------------
library(BDgraph, lib="/storage/groups/epigenereg01/tools/2017/R/3.4")
library(GenomicRanges)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# get snakemake params
# ------------------------------------------------------------------------------
ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]
# iteration/run to be annotated in result table
iteration <- snakemake@params$iteration

print(paste0("Processing file: ", ifile, "."))

# ------------------------------------------------------------------------------
# Perform validation
# ------------------------------------------------------------------------------
# load data
load(ifile)

# the result table
tab <- c()

# check each individual simulated result
temp <- lapply(names(result), function(n) {
  # get validation data
  r <- result[[n]]
  d <- r$data.sim

  # ----------------------------------------------------------------------------
  # Get all fits
  # ----------------------------------------------------------------------------
  ggm_fit <- r$fits$ggm_fit
  ggm_fit_no_priors <- r$fits$ggm_fit_no_priors
  ggm_fit_nor_priors_empty <- r$fits$ggm_fit_no_priors_empty
  genenet_adj <- as(r$fits$genenet_graph,"matrix")
  irn_adj <- as(r$fits$irn_graph,"matrix")

  # ----------------------------------------------------------------------------
  # use the bdgraph internal method to get spec/sens, f1 and MCC. Use the
  # original simulation object containing the ground truth graph
  # ----------------------------------------------------------------------------
  perf <- t(compare(d, ggm_fit, ggm_fit_no_priors, ggm_fit_nor_priors_empty))
  comparisons <- c("true", "ggm_fit", "ggm_fit_no_priors", "ggm_fit_no_priors_empty")
  perf <- as.data.frame(perf)
  rownames(perf) <- paste(n, comparisons, sep="_")

  # same for iRafNet and GeneNet graphs
  perf2 <- t(compare(d, genenet_adj, irn_adj))
  comparisons <- c("true", "genenet_fit", "iRafNet_fit")
  perf2 <- as.data.frame(perf2)
  rownames(perf2) <- paste(n, comparisons, sep="_")
  # only one true graph (we combine in the next line...)
  perf2 <- perf2[!rownames(perf2) %in% "true",]
  perf <- rbind(perf, perf2)

  # remember for easy plotting
  perf$rdegree <- r$rdegree
  perf$snp <- r$snp
  perf$comparison <- comparisons
  perf$iteration <- iteration

  # ----------------------------------------------------------------------------
  # Add to result table
  # ----------------------------------------------------------------------------
  tab <<- rbind(tab, perf)
})

# ------------------------------------------------------------------------------
# write result to output file
write.table(file=ofile, tab, col.names=NA, sep="\t",
            quote=F, row.names = TRUE)

sink()
