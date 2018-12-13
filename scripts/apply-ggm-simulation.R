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

# ------------------------------------------------------------------------------
# Get snakemake input and load data
# ------------------------------------------------------------------------------
fdata <- snakemake@input[[1]]
fout <- snakemake@output[[1]]
threads <- snakemake@threads
nriter <- as.numeric(snakemake@params$nriter)
burnin <- as.numeric(snakemake@params$burnin)
iter <- as.numeric(snakemake@params$iteration)
ntrees <- 1000
npermut <- 5

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
  # number of predictors to be sampled at each node  for the random forest model
  #(we take the suggestion from the package manual for now)
  ntry <- round(sqrt(ncol(d)-1))

  # in case of rbinom simulation, we get adjusted priors
  if(grepl("_rbinom", n)) {
    priors <- sim$priors
  }

  # Start from observed/hidden graph state?
  gstart <- get_gstart_from_priors(priors)

  # ----------------------------------------------------------------------------
  # Get model fits
  # ----------------------------------------------------------------------------
  print("Fitting model using priors.")
  ggm_fit <- bdgraph(d, method = "gcgm",
                     iter = nriter, burnin = burnin,
                     g.prior = priors, g.start = gstart,
                     save.all=T, cores=threads)
  print("Fitting model without priors")
  ggm_fit_no_priors <- bdgraph(d, method = "gcgm",
                               iter = nriter, burnin = burnin,
                               g.start = gstart,
                               save.all=T, cores=threads)
  print("Fitting model without priors using empty start graph.")
  ggm_fit_no_priors_empty <- bdgraph(d, method = "gcgm",
                                           iter = nriter, burnin = burnin,
                                           g.start = "empty",
                                           save.all=T, cores=threads)
  print("Fitting model using iRafNet.")
  irn_out <- iRafNet(d, priors, ntrees, ntry, colnames(d), threads=threads)
  irn_perm_out <- Run_permutation(d, priors,
                                  ntrees, ntry, colnames(d), npermut, threads=threads)
  irn_fit <- iRafNet_network(irn_out, irn_perm_out, TH = 0.05)
  class(irn_fit) <- c(class(irn_fit), "irafnet")

  print("Fitting model using GeneNet.")
  gn_data <- d[,apply(d, 2, function(x) !anyNA(x))]
  pcors <- ggm.estimate.pcor(gn_data)
  genenet_fit <- network.test.edges(pcors, plot = F)
  genenet_fit$node1 <- colnames(gn_data)[genenet_fit$node1]
  genenet_fit$node2 <- colnames(gn_data)[genenet_fit$node2]
  class(genenet_fit) <- c(class(genenet_fit), "genenet")

  # ----------------------------------------------------------------------------
  # Extract graphs and create result lists
  # ----------------------------------------------------------------------------

  # get the result graphs
  ggm_graph <- graph_from_fit(ggm_fit, ranges, annotate=F)
  ggm_graph_no_priors <- graph_from_fit(ggm_fit_no_priors, ranges, annotate=F)
  ggm_graph_no_priors_empty <- graph_from_fit(ggm_fit_no_priors_empty,
                                              ranges, annotate=F)
  irn_graph <- graph_from_fit(irn_fit, ranges, annotate=F)
  genenet_graph <- graph_from_fit(genenet_fit, ranges, annotate=F)

  # new entry in our data collection
  sim$fits <- listN(ggm_fit, ggm_fit_no_priors, ggm_fit_no_priors_empty, irn_fit, genenet_fit,
                    ggm_graph, ggm_graph_no_priors, ggm_graph_no_priors_empty, irn_graph, genenet_graph)
  sim
})

names(result) <- names(simulated_data)

print("Saving results.")
save(file=fout, result, priors, runs)
