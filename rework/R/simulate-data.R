#' Script based on snakemake input/output vars to simulate
#' data for individual sentinels
#' 
#' @author Johann Hawe
#' 

sink(file = snakemake@log[[1]])

#define colors
tropical <- c("darkorange", "dodgerblue", "hotpink", "limegreen", "yellow")
palette(tropical)

# load needed packages
library(graph)
#library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
# cant use central storage for now... remove asap
library(BDgraph)

library(parallel)
# source needed scripts
source("R/priors.R")
source("R/lib.R")

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
create_prior_graphs <- function(priors,
                                rand_steps=seq(0, 1, by=0.1)) {
  
  # now we create randomized prior information
  # for certain degrees
  if(!is.null(rand_steps)) {
    graphs <-list()
    n <- length(rand_steps)
    for(i in 1:n) {
      # sample prior graph
      samp <- sample_prior_graph(priors)
      sentinel <- samp$sentinel
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
#' @param plot.dir Where to save plots to
#' 
#' @return A list of data simulation objects, one for each of the
#' graphs in the input list
#' 
#' @author Johann Hawe
#'
simulate_data <- function(graphs, ggm.data, nodes, plot.dir) {
  
  d <- lapply(graphs, function(g) {
    gr <- g$graph.observed
    s <- g$snp
    
    # create adjacency matrix from our graph
    g_adj <- igraph::as_adj(igraph::igraph.from.graphNEL(gr), sparse = F)
    snp_available <- T
    if(!s %in% colnames(g_adj)) {
      warning("Sentinel not in graph! This should only be the case iff the sentinel has no prior associated.")
      snp_available <- F
    }
    # now we can simulate the data for the entities
    N <- nrow(ggm.data)
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
      gd <- table(round(as.numeric(ggm.data[,s])))
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
ifile <- snakemake@input[[1]]
rwdir <- snakemake@params$random_walk_results
plot.dir <- snakemake@params$plotdir
dir.create(plot.dir)
odir <- snakemake@output[["odir"]]
dir.create(odir)
threads <- snakemake@threads

# load data and utilize lolipop cohort
load(ifile)
# do once before calling the loop
load.gtex.priors()

sentinels <- names(data)

temp <- mclapply(sentinels, function(s) {
  # first snp has only 2 genotype levels, 
  # second one has too many entities for now (takes too long)
  if(s %in% c("rs79755767", "rs60626639")) {
    return(NULL)
  }
  print(paste0("Processing ", s, "..."))
  
  ranges <- data[[s]]$lolipop$ranges
  nodes <- data[[s]]$lolipop$nodes
  ggm.data <- data[[s]]$lolipop$data
  
  priors <- get.link.priors(ranges, nodes)
  
  # create the hidden and observed graphs
  graphs <- create_prior_graphs(priors)
  
  # simulate data for ggm
  simulated_data <- simulate_data(graphs, ggm.data, nodes, plot.dir)
  
  print(paste0("Saving results to ", odir))
  
  # write out the results
  save(file=paste0(odir, s, ".RData"), simulated_data, 
       ranges, nodes, ggm.data)
}, mc.cores=threads)

sink()

