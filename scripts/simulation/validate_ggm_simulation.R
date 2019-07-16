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
library(BDgraph)
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
  bdgraph_fit <- r$fits$bdgraph_fit
  bdgraph_no_priors_fit <- r$fits$bdgraph_no_priors_empty_fit
  glasso_adj <- as(r$fits$glasso, "matrix")
  glasso_no_priors_adj <- as(r$fits$glasso_no_priors, "matrix")
  genenet_adj <- as(r$fits$genenet_graph,"matrix")
  irafnet_adj <- as(r$fits$irafnet_graph,"matrix")

  # ----------------------------------------------------------------------------
  # use the bdgraph internal method to get spec/sens, f1 and MCC. Use the
  # original simulation object containing the ground truth graph
  # ----------------------------------------------------------------------------
  perf <- t(compare(d, bdgraph_fit, bdgraph_no_priors_fit, glasso_adj))
  comparisons <- c("true", "bdgraph", "bdgraph_no_priors", "glasso")
  perf <- as.data.frame(perf)
  rownames(perf) <- paste(n, comparisons, sep="_")

  # same for iRafNet and GeneNet graphs
  perf2 <- t(compare(d, glasso_no_priors_adj, genenet_adj, irafnet_adj))
  comparisons2 <- c("true", "glasso_no_priors", "genenet", "iRafNet")
  perf2 <- as.data.frame(perf2)
  rownames(perf2) <- paste(n, comparisons2, sep="_")
  # only one true graph (we combine in the next line...)
  perf2 <- perf2[!grepl("true", rownames(perf2)),]
  perf <- rbind(perf, perf2)

  # remember for easy plotting
  perf$rdegree <- r$rdegree
  perf$snp <- r$snp
  perf$comparison <- c(comparisons, comparisons2[-1])
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
