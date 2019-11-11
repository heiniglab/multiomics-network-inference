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
library(graph)
library(igraph)
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# get snakemake params
# ------------------------------------------------------------------------------
ffits <- snakemake@input$fits

foutput <- snakemake@output[[1]]

# iteration/run to be annotated in result table
iteration <- snakemake@params$iteration

print(paste0("Processing file: ", ffits, "."))

# ------------------------------------------------------------------------------
# Perform validation
# ------------------------------------------------------------------------------
# load data
load(ffits)

# check each individual simulated result
temp <- lapply(names(result), function(n) {
  # get validation data
  r <- result[[n]]
  d <- r$data.sim

  # ----------------------------------------------------------------------------
  # Get all fits
  # ----------------------------------------------------------------------------
  bdgraph <- r$fits$bdgraph
  bdgraph_no_priors <- r$fits$bdgraph_no_priors
  glasso <- r$fits$glasso
  glasso_no_priors <- r$fits$glasso_no_priors
  genenet <- r$fits$genenet
  irafnet <- r$fits$irafnet
  genie3 <- r$fits$genie3

  gs <- list(bdgraph=bdgraph,
             bdgraph_no_priors=bdgraph_no_priors,
             glasso=glasso,
             glasso_no_priors=glasso_no_priors,
             genenet=genenet,
             irafnet=irafnet,
             genie3=genie3)

  # ----------------------------------------------------------------------------
  # use the bdgraph internal method to get spec/sens, f1 and MCC. Use the
  # original simulation object containing the ground truth graph
  # ----------------------------------------------------------------------------
  perf <- lapply(names(gs), function(g) {
    # we need the adjacency matrix for comparison
    perf <- t(BDgraph::compare(d, as(gs[[g]], "matrix")))
    comparisons <- c("True", g)
    perf <- as.data.frame(perf)
    rownames(perf) <- paste(n, comparisons, sep="_")
    perf <- perf[!grepl("True", rownames(perf)),]

    # annotate density for comparison
    ig <- igraph::igraph.from.graphNEL(gs[[g]])
    dens <- edge_density(ig)
    perf$density_model <- dens
    perf
  })

  perf <- do.call(rbind.data.frame, perf)

  # remember for easy plotting
  perf$rdegree <- r$rdegree
  perf$snp <- r$snp
  perf$iteration <- iteration
  perf$comparison <- names(gs)
  perf$density_true <- edge_density(igraph.from.graphNEL(r$graph.observed))

  perf
})

tab <- do.call(rbind.data.frame, temp)

# ------------------------------------------------------------------------------
# write result to output file
write.table(file=foutput, tab, col.names=NA, sep="\t",
            quote=F, row.names = TRUE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
