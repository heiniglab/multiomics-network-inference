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
  bdgraph_adj <- as(r$fits$bdgraph, "matrix")
  bdgraph_no_priors_adj <- as(r$fits$bdgraph_no_priors, "matrix")
  glasso_adj <- as(r$fits$glasso, "matrix")
  glasso_no_priors_adj <- as(r$fits$glasso_no_priors, "matrix")
  genenet_adj <- as(r$fits$genenet,"matrix")
  irafnet_adj <- as(r$fits$irafnet,"matrix")
  genie3_adj <- as(r$fits$genie3, "matrix")

  gs <- list(bdgraph=bdgraph_adj,
             bdgraph_no_priors=bdgraph_no_priors_adj,
             glasso=glasso_adj,
             glasso_no_priors=glasso_no_priors_adj,
             genenet=genenet_adj,
             irafnet=irafnet_adj,
             genie3=genie3_adj)

  # ----------------------------------------------------------------------------
  # use the bdgraph internal method to get spec/sens, f1 and MCC. Use the
  # original simulation object containing the ground truth graph
  # ----------------------------------------------------------------------------
  perf <- lapply(names(gs), function(g) {
    perf <- t(compare(d, gs[[g]]))
    comparisons <- c("true", g)
    perf <- as.data.frame(perf)
    rownames(perf) <- paste(n, comparisons, sep="_")
    perf <- perf[!grepl("true", rownames(perf)),]
  })

  perf <- do.call(rbind.data.frame, perf)

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
