# ------------------------------------------------------------------------------
#' Script based on snakemake input/output vars to simulate
#' data for individual sentinels
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Prep libraries and source scripts.")
# ------------------------------------------------------------------------------

# load needed packages
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(BDgraph))

# source needed scripts
source("scripts/priors.R")
source("scripts/lib.R")

#define colors
cols <- get_defaultcolors()
palette(cols)

# ------------------------------------------------------------------------------
# Define methods
# ------------------------------------------------------------------------------

#' Creates a graph based on the available priors
#' for a sentinel locus
#'
#' @param priors The (symmetric) prior matrix for the
#' sentinel locus
#' @param rand_steps Vector of randomization rates. Used to randomize
#' prior information with individual difficulty.
#'
#' @return Returns the graph-structure as a graphNEL object
#' based on the priors given as parameters (retaining nodes
#' for which there was no edge)
#'
#' @author Johann Hawe
#'
create_prior_graphs <- function(priors, sentinel,
                                rand_steps=seq(0, 1, by=0.1)) {

  # ----------------------------------------------------------------------------
  # Create distinct randomized graphs
  # ----------------------------------------------------------------------------
  if(!is.null(rand_steps)) {
    graphs <-list()
    n <- length(rand_steps)
    for(i in 1:n) {

      # sample prior graph
      samp <- sample_prior_graph(priors, sentinel)
      gsamp <- samp$sample_graph
      gfull <- samp$full_graph
      ei <- samp$edge_info

      # randomize graph
      rd <- rand_steps[i]

      grand <- rewire_graph(gsamp, rd, ei)
      pe <- get_percent_prioredges(grand$gnew, ei) * 100
      div <- abs((100-pe)-rd*100)

      # report for debug
      #print(paste0("Desired/observed randomization: ", rd*100, "/", 100-pe))

      # check for more than 2 percent divergence in randomization
      if(div >= 2) warning(paste0("Large divergence in randomization: ", div))

      # save graph
      id <- paste0(sentinel, "_rd", rd)
      graphs[[id]] <- list(graph.full=gfull,
                           graph.sampled=gsamp,
                           graph.observed=grand$gnew,
                           rdegree=rd,
                           snp=sentinel)
    }

    # --------------------------------------------------------------------------
    # Now we create a random binomial prior matrix based on graph size
    # --------------------------------------------------------------------------

    # create random binomial pseudo prior for the full graph
    # get a full graph as basis
    gfull <- graphs[[1]]$graph.full
    rd <- "rbinom"
    # fit the binomial on all graph-|E|and get the p
    all_sizes <- unlist(lapply(graphs, function(gr) {
      numEdges(gr$graph.sampled)
    }))

    n <- ncol(priors)
    E <- (n*(n-1))/2

    # use MLE estimator
    prob <- max(sum(all_sizes) / (length(all_sizes) * E), min(priors))

    prbinom <- matrix(prob, ncol=ncol(priors), nrow=nrow(priors))
    colnames(prbinom) <- colnames(priors)
    rownames(prbinom) <- rownames(priors)

    # save graph with prior matrix
    id <- paste0(sentinel, "_rbinom")
    graphs[[id]] <- list(graph.full=gfull,
                         graph.observed=gfull,
                         rdegree=rd,
                         snp=sentinel, priors=prbinom)
    return(graphs)

  } else {
    return(list(graph.full=g, graph.observed=g,
                rd = 0,
                snp=sentinel, priors=priors))
  }
}

#' Method to simulate multivariate data as input for the
#' BDgraph algorithm from a list of given graphNEL objects
#' and the according collected cohort data
#'
#' @param graphs List of graph objects
#' @param ggm.data A data matrix containing data for all entities
#' in the graphs from a specific cohrt (lolipop/kora)
#'
#' @return A list of data simulation objects, one for each of the
#' graphs in the input list
#'
#' @author Johann Hawe
#'
simulate_data <- function(graphs, sentinel, data, nodes) {

  d <- lapply(graphs, function(g) {
    gr <- g$graph.observed
    s <- sentinel

    # create adjacency matrix from our graph
    g_adj <- igraph::as_adj(igraph::igraph.from.graphNEL(gr), sparse = F)
    snp_available <- T
    if(!s %in% colnames(g_adj)) {
      warning("Sentinel not in graph.")
      snp_available <- F
    }
    # now we can simulate the data for the entities
    N <- nrow(data)
    if(length(nodes) != ncol(g_adj)) {
      warning("Number of nodes of collected ggm data and simulated graph differ.")
    }
    p <- ncol(g_adj)

    data.sim <- bdgraph.sim(p, graph=g_adj, N, mean = 0)
    colnames(data.sim$data) <- colnames(g_adj)

    # check whether our SNP data should be available in the simulated
    # data
    if(snp_available) {
      # split the 'gaussian' snp data into 3 groups (0,1,2)
      # for this we utilize allele frequencies of the original
      # genotype data
      gd <- table(round(as.numeric(data[,s])))
      gd <- sort(gd/sum(gd))
      # define needed intervals
      gdi <- c(gd[1], gd[1]+gd[2], 1)

      # transform data into genotypes
      temp <- data.sim$data[,s]
      qs <- (quantile(temp, gdi))
      names(qs) <- names(gd)
      data.sim$data[,s] <- as.numeric(as.character(cut(temp,c(min(temp), qs),
                                                       right = T,
                                                       labels=c("2","1","0"),
                                                       include.lowest = T)))
    }
    g$data.sim <- data.sim
    g
  })
  names(d) <- names(graphs)
  return(d)
}

# ------------------------------------------------------------------------------
# Start processing, get snakemake input
# ------------------------------------------------------------------------------

# input file containing ranges/data of sentinel
fdata <- snakemake@input[["data"]]
franges <- snakemake@input[["ranges"]]
fpriors <- snakemake@input[["priors"]]

# ------------------------------------------------------------------------------
print("Getting snakemake params.")
# ------------------------------------------------------------------------------
fout <- snakemake@output[[1]]
sentinel <- snakemake@params$sentinel
threads <- snakemake@threads
runs <- 1:as.numeric(snakemake@params$runs)

# ------------------------------------------------------------------------------
print("Loading data.")
# ------------------------------------------------------------------------------
data <- readRDS(fdata)
nodes <- colnames(data)
ranges <- readRDS(franges)
priors <- readRDS(fpriors)
# restrict to priors for which we also have data available
priors <- priors[rownames(priors) %in% nodes, colnames(priors) %in% nodes]

# ------------------------------------------------------------------------------
print(paste0("Running ", length(runs), " simulations."))
# ------------------------------------------------------------------------------
simulations <- mclapply(runs, function(x) {
  set.seed(x)

  # create the hidden and observed graphs
  graphs <- create_prior_graphs(priors, sentinel)

  # simulate data for ggm
  simulated_data <- simulate_data(graphs, sentinel, data, nodes)

  simulated_data
}, mc.cores = threads)

names(simulations) <- paste0("run_", runs)

# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
save(file=fout, simulations, priors, ranges, nodes, data,
     runs, fdata, franges, fpriors)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
sink()
sink(type="message")
