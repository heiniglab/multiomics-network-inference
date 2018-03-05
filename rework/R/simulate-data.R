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
#' @param rates Vector of 'FPR/FNR' rates. Used to simulate
#' graphs of individual difficulty. Default is to 
#' 
#' @return Returns the graph-structure as a graphNEL object
#' based on the priors given as parameters (retaining nodes 
#' for which there was no edge)
#' 
#' @author Johann Hawe
#'
create_prior_graphs <- function(priors,
                                rates=seq(0.1, 0.9, by=0.1)) {
  
  # assume pseudo.prior as min-value (should always be
  # the case)
  pseudo.prior <- min(priors)
  
  # the nodes we want to operate on
  nodes <- colnames(priors)
  sentinel <- nodes[grepl("^rs", nodes)]
  
  # create base graph
  g <- graphNEL(nodes, edgemode = "undirected")
  
  # iterate over each combination of nodes
  ee <- expand.grid(nodes, nodes, stringsAsFactors = F)
  colnames(ee) <- c("n1", "n2")
  ee <- (ee[ee$n1 != ee$n2, ])
  # remove duplicated entries
  een <- vector(mode="character", length=nrow(ee))
  for(i in 1:nrow(ee)){
    een[i] <- paste0(sort(c(ee[i,1], ee[i,2])), collapse="_")
  }
  # final collection of all possible edges (undirected graph)
  ee <- ee[!duplicated(een),]
  een <- een[!duplicated(een)]
  rownames(ee) <- een
  
  ee$keep <- F
  ee$prior <- -1
  # iterate over each pair and determine whether to add
  # the pair as an 'true' edge
  set.seed(42)
  for(i in 1:nrow(ee)) {
    # get prior
    p <- priors[ee[i,"n1"], ee[i,"n2"]]
    ee[i,"prior"] <- p
    if(p>pseudo.prior) {
      v <- runif(1)
      if(v<=p) {
        ee[i,"keep"] <- T
      }
    }
  }
  
  # create a full graph for comparison
  keep_full <- ee$prior > pseudo.prior
  gfull <- addEdge(ee[keep_full,1], ee[keep_full,2], g)
  
  # keep only relevent edges and add to graph
  to_add <- ee[ee$keep,,drop=F]
  if(nrow(to_add)>0) {
    # make sure we add the SNP
    # -> if its not in the list of edges to keep,
    # we simply switch one random connection with it
    # if it doesn't have a prior, though, we just ignore it
    if(!(sentinel %in% c(to_add[,"n1"], to_add[,"n2"]))) {
      temp <- priors[,sentinel]
      temp <- temp[temp>pseudo.prior]
      if(!length(temp) == 0) {
        # select random name
        temp2 <- sample(names(temp),1)
        to_switch <- sample(1:nrow(to_add),1)
        old <- to_add[to_switch,]
        old["keep"] <- F
        
        to_add[to_switch,] <- c(sentinel, temp2, T, priors[sentinel,temp2])
        rownames(to_add)[to_switch] <- paste0(sort(c(sentinel, temp2)), collapse="_")
        ee[ee$n1 == old$n1 & ee$n2 == old$n2,] <- old
      }
    }
  } else {
    # none of the edges made it
    # for now add at least the edges which got any prior so that we do not 
    # get an empty graph here; add random SNP edge to have the SNP available
    ee[ee$prior>pseudo.prior,"keep"] <- T
    ee[which(grepl("^rs", ee$n1) | grepl("^rs", ee$n2))[1], "keep"] <- T
    to_add <- ee[ee$keep,,drop=F]
  }
  
  g <- addEdge(to_add[,1], to_add[,2], g)
  
  # now we add the 'false prior' edges by rewiring
  if(!is.null(rates)) {
    graphs <-list()
    n <- length(rates)
    for(i in 1:n) {
      # the proportion of prior edges to keep
      prior_prop <- rates[i]
      needed_prior_edges <- round(prior_prop * numEdges(g))
      g_rewired <- rewire_graph(g, 1)
      edges_in_priors <- numEdges(graph::intersection(g_rewired, gfull))
      counter <- 1
      # we try at most 100 rewirings, otherwise we use the one closest to our ratio
      max <- 100000
      g_temp <- g_rewired
      print("Started rewiring, this might take a while.")
      while(needed_prior_edges != edges_in_priors)
      {
        to_rewire <- edges_in_priors - needed_prior_edges + 1
        g_rewired <<- rewire_graph(g_rewired, to_rewire)
        n <- numEdges(graph::intersection(g_rewired, gfull))
        if(abs(needed_prior_edges - n) < abs(needed_prior_edges - edges_in_priors)) {
          g_temp <- g_rewired
          edges_in_priors <- n
        }
        counter <- counter + 1
        if(counter %% 100 == 0) {
          print(paste0(counter, "  ", n, "  ", needed_prior_edges))
        }
        if(counter == max) {
          break
        }
      }
      print(paste0("Done.\nTried ", counter, " rewirings."))
        
      # filter zero degree nodes from the graphs
      g <- remove_unconnected_nodes(g)
      g_rewired <- remove_unconnected_nodes(g_rewired)
      
      # drop nodes without edges
      # save graph
      id <- paste0(sentinel, "_fpr", drop_prop)
      graphs[[id]] <- list(graph.hidden=g, 
                        graph.observed=gr,
                        prior_rate=fpr, fnr=fnr,
                        snp=sentinel, priors=priors)
    }
    return(graphs)
  } else {
    return(list(graph.hidden=g, graph.observed=g, 
                fpr=0, fnr=0,
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
    # # compare7plot the two allele frequencies
    # gd2 <- table(data.sim$data[,s])
    # gd2 <- sort(gd2/sum(gd2))
    # 
    # pdf(paste0(plot.dir, "/", s, "-genotype-frequencies.pdf"))
    # barplot(gd, col=1, main=paste0(s, " genotype frequencies"))
    # cs <- col2rgb(2)[,1]/255
    # cs <- rgb(cs[1], cs[2], cs[3],
    #           0.5)
    # barplot(gd2, add = T, col=cs)
    # legend("topleft", legend = c("original", "transformed"), fill=c(1,cs))
    # dev.off()
    g$data.sim <- data.sim
    g
  })
  names(d) <- names(graphs)
  return(d)
}

rewire_graph <- function(g, prop) {
  
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
  
  # we only need the plotgraph from the 
  # random walk result
  # currently not available!
  #env <- new.env()
  #rwfile <- paste0(rwdir, s, ".RData")
  #if(!file.exists(rwfile)) {
  #  return(NULL)
  #}
  #load(rwfile, env)
  #rw_graph <- with(env, plot.graph)
  #rw_graph_nodes <- nodes(rw_graph)
  
  ranges <- data[[s]]$lolipop$ranges
  nodes <- data[[s]]$lolipop$nodes
  ggm.data <- data[[s]]$lolipop$data

  # we add some random edges based on 
  # the nodes in our data to the ground truth
  # first we add all remaining nodes
 # toadd <- setdiff(nodes, rw_graph_nodes)
#  gr <- graph::addNode(toadd, rw_graph)
  
  # now add some random edges between the newly added nodes
  # and the nodes which were already in the graph
  set.seed(42)
  # the max number of edges to be randomly added to the graph
  #r <- numEdges(rw_graph)
  #gr <- addEdge(sample(rw_graph_nodes, r, replace = T), 
  #              sample(toadd, r, replace = T), gr)
  
  # annotate our graph with the appropriate node- and edgeData
  #gr <- annotate.graph(gr, ranges)
  #gr <- list(graph.hidden=rw_graph, 
  #     graph.observed=gr,
  #     fpr=-1, fnr=-1)
  
  # get priros for prior based graphs
  priors <- get.link.priors(ranges, nodes)
  
  # create the hidden and observed graphs
  graphs <- create_prior_graphs(priors)
  gn <- names(graphs)
 #graphs <- append(graphs, gr)
 #names(graphs) <- c(gn, "random_walk")
  
  # simulate data for ggm
  simulated_data <- simulate_data(graphs, ggm.data, nodes, plot.dir)
  
  print(paste0("Saving results to ", odir))
  
  # write out the results
  save(file=paste0(odir, s, ".RData"), simulated_data, 
       ranges, nodes, ggm.data)
}, mc.cores=threads)

sink()

