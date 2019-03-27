#' -----------------------------------------------------------------------------
#' Apply the network inference using different models on simulated data
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Mar 26 10:45:32 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(GenomicRanges))
library(pheatmap)
library(doParallel)
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(graph))

source("scripts/lib.R")
source("scripts/reg_net.R")

# ------------------------------------------------------------------------------
print("Get snakemake input and load data.")
# ------------------------------------------------------------------------------
# inputs
fdata <- snakemake@input$data
fppi_db <- snakemake@input$ppi_db
fcpg_context <- snakemake@input$cpg_context

# outputs
fout <- snakemake@output[[1]]

# params
threads <- snakemake@threads
sim_iter <- as.numeric(snakemake@params$iteration)

# register clusters for multi-threading (used in iRafNet)
cl <- makeCluster(threads)
registerDoParallel(cl)

# contains: simulations, ranges, priors, nodes, data, run
load(fdata)
ppi_db <- readRDS(fppi_db)

# ------------------------------------------------------------------------------
# we generated several graphs, which we all pcalculate models for now
# ------------------------------------------------------------------------------

# apply over the different runs/iterations
simulated_data <- simulations[[sim_iter]]

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

  # ------------------------------------------------------------------------------
  print("Infer regulatory networks.")
  # ------------------------------------------------------------------------------
  print("Fitting model using glasso.")
  gl <- glasso_screen(d, priors, threads, ranges, ppi_db, fcpg_context)
  glasso <- gl$glasso
  glasso_no_priors <- gl$glasso_no_priors
  glasso_all <- gl$glasso_all
  print("Done with glasso screen.")

  print("Fitting model using priors.")
  bdgraph <- reg_net(d, priors, "bdgraph", threads=threads)

  print("Fitting model without priors using empty start graph.")
  bdgraph_no_priors_empty <- reg_net(d, NULL, "bdgraph",
                                     use_gstart = F, threads=threads)

  print("Fitting model using iRafNet.")
  irafnet <- reg_net(d, priors, "irafnet", threads=threads)

  print("Fitting model using GeneNet.")
  genenet <- reg_net(d, priors, "genenet", threads=threads)

  # ------------------------------------------------------------------------------
  print("Add custom annotations for the graphs.")
  # ------------------------------------------------------------------------------
  bdgraph$graph <- annotate.graph(bdgraph$graph, ranges, ppi_db, fcpg_context)
  bdgraph_no_priors_empty$graph <- annotate.graph(bdgraph_no_priors_empty$graph,
                                                  ranges, ppi_db, fcpg_context)
  irafnet$graph <- annotate.graph(irafnet$graph, ranges, ppi_db, fcpg_context)
  genenet$graph <- annotate.graph(genenet$graph, ranges, ppi_db, fcpg_context)

  # create result list
  result <- list(bdgraph_fit = bdgraph$fit,
                 bdgraph_no_priors_empty_fit = bdgraph_no_priors_empty$fit,
                 irafnet_fit = irafnet$fit,
                 genenet_fit = genenet$fit,
                 bdgraph = bdgraph$graph,
                 bdgraph_no_priors_empty = bdgraph_no_priors_empty$graph,
                 irafnet_graph = irafnet$graph,
                 genenet_graph = genenet$graph,
		 glasso_fit = glasso$fit,
		 glasso = glasso$graph,
		 glasso_no_priors_fit = glasso_no_priors$fit,
		 glasso_no_priors = glasso_no_priors$graph,
		 glasso_all = glasso_all)

  sim$fits <- result

  sim
})

names(result) <- names(simulated_data)

print("Saving results.")
save(file=fout, result, priors, runs)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
