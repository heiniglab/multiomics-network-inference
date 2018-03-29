#' Script based on snakemake input/output vars to simulate
#' data for individual sentinels
#' 
#' @author Johann Hawe
#' 

#define colors
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)

# load needed packages
library(graph)
library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4/")

# source needed scripts
source("scripts/priors.R")
source("scripts/lib.R")

#' Creates a graph based on the available priors
#' for a sentinel locus
#'
#' @param priors The (symmetric) prior matrix for the 
#' sentinel locus
#' @param rates Vector of randomization rates. Used to randomize
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
  
  # now we create randomized prior information
  # for certain degrees
  if(!is.null(rand_steps)) {
    graphs <-list()
    n <- length(rand_steps)
    for(i in 1:n) {
      # sample prior graph
      samp <- sample_prior_graph(priors, sentinel)
      g <- samp$sample_graph
      gfull <- samp$full_graph
      
      # the proportion of prior edges to keep
      # 'randomization degree'
      rd <- rand_steps[i]
        
      # get randomized priors
      prand <- randomize_priors(priors, rd)
      
      # filter zero degree nodes from the graphs
      #g <- remove_unconnected_nodes(g)
      #gfull <- remove_unconnected_nodes(gfull)
      
      # save graph
      id <- paste0(sentinel, "_rd", rd)
      graphs[[id]] <- list(graph.full=gfull, 
                        graph.observed=g,
                        rdegree=rd,
                        snp=sentinel, priors=prand)
    }
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

# input file containing ranges/data of sentinel
fdata <- snakemake@input[["data"]]
franges <- snakemake@input[["ranges"]]
fpriors <- snakemake@input[["priors"]]
fout <- snakemake@output[[1]]
sentinel <- snakemake@params$sentinel
threads <- snakemake@threads

# load data
data <- readRDS(fdata)
ranges <- readRDS(franges)
priors <- readRDS(fpriors)

# do the simulation 100 times
simulations <- mclapply(1:100, function(x) {
  print(paste0("Number of runs performed:", x))
  set.seed(x)
  # create the hidden and observed graphs
  graphs <- create_prior_graphs(priors, sentinel)
  nodes <- colnames(data)
  
  # simulate data for ggm
  simulated_data <- simulate_data(graphs, sentinel, data, nodes)
  
  simulated_data
}, mc.cores = threads)

names(simulations) <- paste0("run_", 1:100)
# write out the results
save(file=fout, simulations, ranges, nodes, data)