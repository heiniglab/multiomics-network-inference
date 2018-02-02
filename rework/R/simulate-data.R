#' Script based on snakemake input/output vars to simulate
#' data for individual sentinels
#' 
#' @author Johann Hawe
#' 

sink(file = snakemake@log[[1]])

# load needed packages
library(graph)
library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4/")
# source needed scripts
source("R/priors.R")
source("R/lib.R")

# input file containing ranges/data of sentinel
ifile <- snakemake@input[[1]]
rwdir <- snakemake@params$random_walk_results
odir <- snakemake@output[["odir"]]

# load data and utilize lolipop cohort
load(ifile)

for(s in names(data)) {
  print(paste0("Processing ", s, "..."))
  
  # we only need the plotgraph from the 
  # random walk result
  env <- new.env()
  rwfile <- paste0(rwdir, s, ".RData")
  if(!file.exists(rwfile)) {
    next
  }
  
  load(rwfile, env)
  rw_graph <- with(env, plot.graph)
  rw_graph_nodes <- nodes(rw_graph)
  
  ranges <- data[[s]]$lolipop$ranges
  nodes <- data[[s]]$lolipop$nodes
  ggm.data <- data[[s]]$lolipop$data

  # we add some random edges based on 
  # the nodes in our data to the ground truth
  # first we add all remaining nodes
  toadd <- setdiff(nodes, rw_graph_nodes)
  gr <- addNode(toadd, rw_graph)
  
  # now add some random edges between the newly added nodes
  # and the nodes which were already in the graph
  set.seed(42)
  # the max number of edges to be randomly added to the graph
  r <- ceiling(length(rw_graph_nodes) / 2)
  r <- min(r, length(toadd))
  gr <- addEdge(sample(rw_graph_nodes, r), sample(toadd, r), gr)
  
  # annotate our graph with the appropriate node- and edgeData
  gr <- annotate.graph(gr, ranges)

  # create adjacency matrix from our graph
  gr_adj <- igraph::as_adj(igraph::igraph.from.graphNEL(gr), sparse = F)
  # create a 'well-informed' covariance and precision matrix
  # TODO we likely not need it, but as it is now it is not working
  # at all!
#  sigma <- priors
#  # here we add some noise on the priors
#  len <- dim(sigma)[1] * dim(sigma)[2]
#  noise <- matrix(rnorm(len, 0, 1), dim(sigma)[1])
#  sigma <- sigma + noise
#  # diag should be all 1s, since this is cov(x,x)
#  diag(sigma) <- 1
#  K <- solve(sigma)
#  # define true graph
#  # for now just use the start graph (i.e. where there was any
#  # prior information available)
#  G <- get.g.start.from.priors(priors)
  
  # now we can simulate the data for the entities
  N <- nrow(ggm.data)
  if(length(nodes) != ncol(gr_adj)) {
    warning("Number of nodes of collected ggm data and simulated graph differ.")
  }
  p <- ncol(gr_adj)

  data.sim <- bdgraph.sim(p, graph=gr_adj, N, mean = 0)
  colnames(data.sim$data) <- colnames(gr_adj)
  # split the 'gaussian' snp data into 3 groups (0,1,2)
  # TODO do this more sophistcally! this is but a quick hack...
  data.sim$data[, s] <- as.numeric(cut(data.sim$data[, s], 
                                       breaks = 3, 
                                       labels=c(0,1,2)))
  
  print(paste0("Saving results to ", odir))
  
  # write out the results
  save(file=paste0(odir, s, ".RData"), data.sim, N, p, gr, 
       gr_adj, ranges, nodes, ggm.data)
}

sink()

