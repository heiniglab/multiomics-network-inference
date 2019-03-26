# ------------------------------------------------------------------------------
# load needed libraries and source scripts
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(BDgraph,
                                       lib="/storage/groups/epigenereg01/tools/2017/R/3.4/"))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(iRafNet))
suppressPackageStartupMessages(library(GeneNet))

source("scripts/lib.R")
source("scripts/bdgraph-supplement.R")
source("scripts/reg_net.R")

# ------------------------------------------------------------------------------
# Get snakemake input and load data
# ------------------------------------------------------------------------------
fdata <- snakemake@input[[1]]
fout <- snakemake@output[[1]]
threads <- snakemake@threads

# number of simulation iterations performed
iter <- as.numeric(snakemake@params$iteration)

# currently not in used!
#nriter <- as.numeric(snakemake@params$nriter)
#burnin <- as.numeric(snakemake@params$burnin)

#ntrees <- 1000
#npermut <- 5

# contains: simulations, ranges, priors, nodes, data, run
load(fdata)

# ------------------------------------------------------------------------------
# we generated several graphs, which we all pcalculate models for now
# ------------------------------------------------------------------------------

# apply over the different runs/iterations
simulated_data <- simulations[[iter]]

# for this run, apply over all simulated graphs (different randomization
# degrees)
result <- lapply(names(simulated_data), function(n) {
  # ----------------------------------------------------------------------------
  # Get data
  # ----------------------------------------------------------------------------
  sim <- simulated_data[[n]]
  # sentinel name and simulated data
  s <- sim$snp
  d <- sim$data.sim$data

  # in case of rbinom simulation, we get adjusted priors
  if(grepl("_rbinom", n)) {
    priors <- sim$priors
  }

  # obtain a start graph
  gstart <- get_gstart_from_priors(priors)

  # ----------------------------------------------------------------------------
  # Get model fits
  # ----------------------------------------------------------------------------
  print("Fitting model using priors.")
  bdgraph <- reg_net(d, priors, "bdgraph", threads=threads)

  print("Fitting model without priors without start graph")
  bdgraph_empty <- reg_net(d, priors, "bdgraph",
                               use_gstart = F, threads=threads)

  print("Fitting model without priors with start graph")
  bdgraph_no_priors <- reg_net(d, NULL, "bdgraph",
                               gstart = gstart, threads=threads)

  print("Fitting model without priors using empty start graph.")
  bdgraph_no_priors_empty <- reg_net(d, NULL, "bdgraph",
                                     use_gstart = F, threads=threads)

  print("Fitting model using iRafNet.")
  irafnet <- reg_net(d, priors, "irafnet", threads=threads)

  print("Fitting model using GeneNet.")
  genenet <- reg_net(d, priors, "genenet", threads=threads)

  # create result list
  result <- list(ggm_fit = bdgraph$fit,
                 ggm_fit_empty = bdgraph_empty$fit,
                 ggm_fit_no_priors = bdgraph_no_priors$fit,
                 ggm_fit_no_priors_empty = bdgraph_no_priors_empty$fit,
                 irn_fit = irafnet$fit,
                 genenet_fit = genenet$fit,
                 ggm_graph = bdgraph$graph,
                 ggm_graph_empty = bdgraph_empty$graph,
                 ggm_graph_no_priors = bdgraph_no_priors$graph,
                 ggm_graph_no_priors_empty = bdgraph_no_priors_empty$graph,
                 irn_graph = irafnet$graph,
                 genenet_graph = genenet$graph)

  sim$fits <- result

  sim
})

names(result) <- names(simulated_data)

print("Saving results.")
save(file=fout, result, priors, runs)
