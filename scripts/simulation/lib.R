#' -----------------------------------------------------------------------------
#' Library of methods for simulation study
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Mar 26 10:04:26 2020
#' -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Method to quickly rewire a given graph
#'
#' @param g The graphNEL object to be rewired
#' @param p The probability with which an edge gets rewired
#' @param ei Edge info object
#'
#' @return The rewired graph object
#'
# ------------------------------------------------------------------------------
rewire_graph <- function(g, p, ei) {
  if(!"edge_info" %in% class(ei)){
    stop("No edge info object provided")
  }
  
  # ----------------------------------------------------------------------------
  # get node information (degrees, non_prior nodes)
  degs <- cbind.data.frame(sort(graph::degree(g), decreasing = T))
  colnames(degs) <- c("degree")
  degs$node <- rownames(degs)
  
  # number of edges to rewire
  N <- round(p * numEdges(g))
  # do nothing...
  if(N == 0) {
    return(list(gnew=g, to_switch=NULL))
  }
  
  # get the 0-degree nodes which do not have any prior associated
  ei_priors <- ei[ei$keep_full,,drop=F]
  nodes <- nodes(g)
  nodes_no_priors <- nodes[!(nodes %in% ei_priors$n1) &
                             !(nodes %in% ei_priors$n2)]
  
  # ----------------------------------------------------------------------------
  # select nodes to switch such that sum of node degree == N
  idxs <- which(degs$degree > 0)
  # maximum Nmax non-prior nodes are available for switching
  Nmax <- length(nodes_no_priors)
  lo <- N-1
  up <- N+1
  # the best node combination (least distance to N)
  best <- NULL
  best_diff <- 0.5
  # we save any found combination of nodes
  # in case we dont have enough prior nodes
  last <- NULL
  add_rewiring <- F
  
  cnt <- 0
  while(best_diff > 0 & cnt < 500) {
    cnt <- cnt + 1
    res <- degreesum_in_range(degs, sample(idxs),
                              lo, up,
                              g)
    
    if(!is.null(res) & !is.null(res$nodes)) {
      last <- res$nodes
      s <- sum(last$degree)
      a <- res$adjustment
      # get 'best' result -> closest to desired number of change
      d <- abs(s-a-N)
      n <- nrow(last)
      if(d < best_diff & n <= Nmax) {
        best <- last
        best_diff <- d
      }
    }
  }
  # check whether we need to do additional work
  if(is.null(best)) {
    best <- last
    # this might be a significant change! was (falsely) "F" before
    add_rewiring <- T
  }
  
  # ----------------------------------------------------------------------------
  # Perform switching
  Ns <- nrow(best)
  
  if(length(nodes_no_priors) > 0) {
    # define switchings
    Nnp <- length(nodes_no_priors)
    
    to_switch <- cbind.data.frame(n_prior=best$node,
                                  n_prior_degree=best$degree,
                                  n_noprior=sample(nodes_no_priors,
                                                   size = Ns,
                                                   replace = Ns>Nnp),
                                  stringsAsFactors=F)
    
    gnew <- switch_nodes(g, to_switch)
    
    # check whether we need to do additional rewiring
    div <- get_prior_divergence(gnew, ei, p)
    if(add_rewiring & div > 2) {
      cnt <- 0
      ignew <- igraph.from.graphNEL(gnew, weight=F)
      
      while(div > 2 & cnt < 500) {
        cnt <- cnt + 1
        ignew <- rewire(ignew, igraph::keeping_degseq(niter=1))
        div <- get_prior_divergence(ignew, ei, p)
      }
      
      gnew <- igraph.to.graphNEL(ignew)
    }
    
  } else {
    stop("Not implemented yet: There where no nodes without priors")
  }
  
  return(listN(gnew, to_switch))
}

#' -----------------------------------------------------------------------------
#' Gets a percent value [0,100] of how much the observed prior edges in a graph
#' diverge from the desired fraction of prior edges.
#'
#' @param g The graph to be checked
#' @param ei The corresponding edge_info object
#' @param p The percentage of desired randomized graph edges
#'
#' @author Johann Hawe
#' -----------------------------------------------------------------------------
get_prior_divergence <- function(g, ei, p) {
  pe <- get_percent_prioredges(g, ei) * 100
  div <- abs((100-pe)-p*100)
  return(div)
}

#' -----------------------------------------------------------------------------
#' Creates a graph with switched nodes as compared to a base graph
#'
#' @param g The graphNEL object in which to switch the nodes
#' @param to_switch The dataframe of nodes to switch (prior vs noprior nodes)
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
switch_nodes <- function(g, to_switch) {
  # get the full edge matrix and switch nodes
  n <- nodes(g)
  em <- t(edgeMatrix(g))
  em <- cbind.data.frame(n1=n[em[,1]],
                         n2=n[em[,2]],
                         stringsAsFactors=F)
  
  for(i in 1:nrow(to_switch)) {
    prior_node <- to_switch[i,"n_prior"]
    noprior_node <- to_switch[i,"n_noprior"]
    # replace all prior node occurrences with the nonprior node
    em[em$n1 == prior_node,"n1"] <- noprior_node
    em[em$n2 == prior_node,"n2"] <- noprior_node
  }
  
  g_new <- graphNEL(n)
  g_new <- addEdge(em[,1], em[,2], g_new)
  return(g_new)
}

#' -----------------------------------------------------------------------------
#' Gets a set of nodes from a dataframe with annotated degree,
#' where the sum of all degrees is within a certain window.
#'
#' @param df The dataframe containing the nodes and their degrees (degree column)
#' @param idxs Indicies of the nodes in the df to be used
#' @param lo Lower bound for the sum
#' @param up Upper bound for the sum
#' @param g The graph from which the nodes originate. If given, it is ensured
#' that nodes in the final collection of nodes are not neighbours within the graph
#'
#' @return Set of nodes for which the sum of degrees is between lo and up
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
degreesum_in_range <- function(df, idxs, lo, up, g) {
  
  # get the  center value for exact matching
  total <- (lo+up)/2
  
  # check special cases
  if(total == 0) {
    return(NULL)
  }
  
  rand <- c()
  # counter for total running sum of degrees we added
  running_sum <- 0
  # counter of substracted information, i.e. whenever we find a new node has
  # neighbours within our collection, we reduce the total amount of degrees found
  # by the respective number of neighbours. We need this information for
  # the evaluation of the found set of nodes
  total_subst <- 0
  
  for(i in idxs) {
    pi <- df[i,,drop=F]
    pi <- cbind(pi, idx=i)
    v <- pi[,"degree"]
    s <- running_sum + v
    subst <- 0
    # neighbour already in collection? handle that
    # adjust total sum of degrees we switch when using this one
    if(!is.null(nrow(rand))) {
      nn <- rand[,"node"]
      # get all neighbours of current node to check how many are already in our
      # collection
      cn <- pi[,"node"]
      neigh <- unlist(graph::adj(g, cn))
      subst <- sum(nn %in% neigh)
      s <- s - subst
    }
    # did not reach lower bound yet
    if(s<lo) {
      running_sum <- s
      total_subst <- total_subst + subst
      rand <- rbind(rand,pi)
      next
    }
    # above lower bound
    if(s>=lo) {
      # match exactly? then done
      if(s==total) {
        running_sum <- s
        total_subst <- total_subst + subst
        rand <- rbind(rand,pi)
        break
      }
      # below upper bound
      if(s<=up) {
        rand <- rbind(rand,pi)
        running_sum <- s
        # first above exact? then done
        if(s>total) {
          running_sum <- s
          total_subst <- total_subst + subst
          break
        }
      } else {
        # exceeded with current value, stop
        running_sum <- s
        total_subst <- total_subst + subst
        break
      }
    }
  }
  
  return(list(nodes=rand, adjustment=total_subst))
}

#' -----------------------------------------------------------------------------
#' Helper method calculating the percentage of edges in the given graph object
#' which have a prior according to the given prior matrix
#'
#' @param g The graph for which to check the edges, igraph or graphNEL object
#' @param ei An edge_info object containing all prior information
#'
#' @return Percentage of edges in the graph which have a prior
#'
#' -----------------------------------------------------------------------------
get_percent_prioredges <- function(g, ei) {
  
  library(graph)
  
  # get prior edges
  # all 'keep_full' edges have a prior > pseudo.prior
  prior_edges <- ei[ei$keep_full,,drop=F]
  ei_nodes <- unique(c(ei$n1, ei$n2))
  
  # ----------------------------------------------------------------------------
  # get the edgematrix with nodes
  if(class(g) == "igraph") {
    n <- names(igraph::V(g))
  } else {
    n <- graph::nodes(g)
  }
  # sanity check
  if(!all(n %in% ei_nodes)) {
    warning("Not all nodes present in edge_info.")
    return(NULL)
  }
  
  # get the edge information from the graph
  if(class(g) == "igraph") {
    em <- igraph::as_edgelist(g)
  } else {
    em <- t(edgeMatrix(g))
  }
  if(nrow(em) > 0) {
    count <- 0
    if(class(g) == "graphNEL") {
      em <- cbind.data.frame(n[em[, 1]],
                             n[em[, 2]],
                             stringsAsFactors = F)
    }
    
    # --------------------------------------------------------------------------
    # check for the number of edges of g being in the prior edges
    for(i in 1:nrow(em)) {
      gn1 <- em[i,1]
      gn2 <- em[i,2]
      
      # check in edge_info prior object
      temp <- any((prior_edges$n1 %in% gn1 & prior_edges$n2 %in% gn2) |
                    (prior_edges$n1 %in% gn2 & prior_edges$n2 %in% gn1))
      
      if(temp) count <- count + 1
    }
    return(count/nrow(em))
  } else {
    warning("No edges in graph.")
    return(NULL)
  }
}

#' -----------------------------------------------------------------------------
#' Samples a graph from the given prior matrix
#' Uses the colnames of the (symmetric) given prior
#' matrix to determine graph nodes
#'
#' @param priors The prior matrix
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
sample_prior_graph <- function(priors, sentinel) {
  
  library(reshape2)
  
  # assume pseudo.prior as min-value (should always be
  # the case)
  pseudo.prior <- min(priors)
  
  # the nodes we want to operate on
  nodes <- colnames(priors)
  
  # create base graph
  g <- graphNEL(nodes, edgemode = "undirected")
  
  # iterate over each combination of nodes
  temp <- priors
  temp[upper.tri(temp, T)] <- NA
  ee <- melt(temp, na.rm=T, stringsAsFactors=F)
  colnames(ee) <- c("n1", "n2", "prior")
  ee$n1 <- as.character(ee$n1)
  ee$n2 <- as.character(ee$n2)
  
  # create names
  rownames(ee) <- paste(ee$n1, ee$n2, sep="_")
  
  # flag for sampling
  ee$keep <- F
  
  # iterate over each pair and determine whether to add
  # the pair as a 'true' edge
  for(i in 1:nrow(ee)) {
    # get prior
    p <- ee[i,"prior"]
    if(p>pseudo.prior) {
      v <- runif(1)
      if(v<=p) {
        ee[i,"keep"] <- T
      }
    }
  }
  
  # create a full graph for comparison
  ee$keep_full <- ee$prior > pseudo.prior
  gfull <- addEdge(ee[ee$keep_full,1], ee[ee$keep_full,2], g)
  
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
        
        to_add[to_switch,] <- c(sentinel, temp2, priors[sentinel,temp2], T, T)
        rownames(to_add)[to_switch] <- paste0(sort(c(sentinel, temp2)),
                                              collapse="_")
        ee[ee$n1 == old$n1 & ee$n2 == old$n2,] <- old
      }
    }
  } else {
    # none of the edges made it
    # for now add at least the edges which got any prior so that we do not
    # get an empty graph here; add random SNP edge to have the SNP available
    warning("Creating full graph as observed graph.")
    ee[ee$prior>pseudo.prior,"keep"] <- T
    ee[which(grepl("^rs", ee$n1) | grepl("^rs", ee$n2))[1], "keep"] <- T
    to_add <- ee[ee$keep,,drop=F]
  }
  
  g <- addEdge(to_add[,1], to_add[,2], g)
  class(ee) <- c(class(ee), "edge_info")
  return(list(sample_graph=g, full_graph=gfull,
              edge_info=ee, sentinel=sentinel))
}

#' -----------------------------------------------------------------------------
#' Creates a graph based on the available priors for a sentinel locus.
#' 
#' Strategy is as follows: For each randomization step (rand_steps), sample
#' a prior graph (uniformly from priors).
#'
#' @param priors The (symmetric) prior matrix for the
#' sentinel locus
#' @param sentinel ID of the sentinel for the blueprint hotspot
#' @param rand_steps Vector of randomization rates. Used to randomize
#' prior information with individual difficulty.
#'
#' @return Returns the graph-structure as a graphNEL object
#' based on the priors given as parameters (retaining nodes
#' for which there was no edge)
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
create_prior_graphs <- function(priors, sentinel,
                                rand_steps=seq(0, 1, by=0.1),
                                threads = 1) {
  require(parallel)
  
  # ----------------------------------------------------------------------------
  # Create distinct randomized graphs
  if(!is.null(rand_steps)) {
    n <- 1:length(rand_steps)
    
    graphs <- mclapply(n, function(i) {
     # sample prior graph
      samp <- sample_prior_graph(priors, sentinel)
      gsamp <- samp$sample_graph
      gfull <- samp$full_graph
      ei <- samp$edge_info
      
      # randomize graph
      rd <- rand_steps[i]
      
      grand <- rewire_graph(gsamp, rd, ei)
      pe <- get_percent_prioredges(grand$gnew, ei) * 100
      div <- abs((100 - pe) - rd * 100)
     
      # check for more than 2 percent divergence in randomization
      if (div >= 2)
        warning(paste0("Large divergence in randomization: ", div))
      
      # save graph
      g <- list(
        graph.full = gfull,
        graph.sampled = gsamp,
        graph.observed = grand$gnew,
        rdegree = rd,
        snp = sentinel
      )
      g
    }, mc.cores = threads) # ---------------------------------------------------
    names(graphs) <- paste0(sentinel, "_rd", rand_steps)
    
    # --------------------------------------------------------------------------
    # Now we create a random binomial prior matrix based on graph size
    
    # create random binomial pseudo prior for the full graph
    # get a full graph as basis
    gfull <- graphs[[1]]$graph.full
    rd <- "rbinom"
    
    # fit the binomial on all graph-|E|and get the p
    all_sizes <- unlist(lapply(graphs, function(gr) {
      numEdges(gr$graph.sampled)
    }))
    
    n <- ncol(priors)
    E <- (n * (n - 1)) / 2
    
    # use MLE estimator
    prob <-
      max(sum(all_sizes) / (length(all_sizes) * E), min(priors))
    
    prbinom <- matrix(prob, ncol = ncol(priors), nrow = nrow(priors))
    colnames(prbinom) <- colnames(priors)
    rownames(prbinom) <- rownames(priors)
    
    # save graph with prior matrix
    id <- paste0(sentinel, "_rbinom")
    graphs[[id]] <- list(
      graph.full = gfull,
      graph.observed = gfull,
      rdegree = rd,
      snp = sentinel,
      priors = prbinom
    )
    
    return(graphs)
    
  } else {
    return(list(graph.full=g, 
                graph.observed=g,
                rd = 0,
                snp=sentinel, 
                priors=priors))
  }
}

#' -----------------------------------------------------------------------------
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
#' -----------------------------------------------------------------------------
simulate_data <- function(graphs, sentinel, data, nodes, threads = 1) {
  require(parallel)
  
  d <- mclapply(graphs, function(g) {
    gr <- g$graph.observed
    s <- sentinel
    
    # create adjacency matrix from our graph
    g_adj <- as(gr, "matrix")
    
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
      if(length(gd) == 3) {
        # define needed intervals
        gdi <- c(gd[1], gd[1] + gd[2], 1)
        labels <- c("2", "1", "0")
      } else {
        # define needed intervals
        gdi <- c(gd[1], gd[1] + gd[2])
        labels <- names(gd)
      }
      
      # transform data into genotypes
      temp <- data.sim$data[,s]
      qs <- (quantile(temp, gdi))
      names(qs) <- names(gd)
      data.sim$data[,s] <- as.numeric(as.character(cut(temp, c(min(temp), qs),
                                                       right = T,
                                                       labels=labels,
                                                       include.lowest = T)))
    }
    g$data.sim <- data.sim
    g
  }, mc.cores = threads) # -----------------------------------------------------
  names(d) <- names(graphs)
  return(d)
}

# ------------------------------------------------------------------------------
#' Generate the validation table for all simul results contained in the provided
#' data object.
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' 
# ------------------------------------------------------------------------------
get_validation_table <- function(result, iteration) {
  
  require(dplyr)
  require(igraph)
  
  # check each individual simulated result
  tab <- lapply(names(result), function(n) {
    # get validation data
    r <- result[[n]]
    d <- r$data.sim
    
    # --------------------------------------------------------------------------
    # Get all fitted graphs, remove the corresponding 'fit' objects
    gs <- r$fits[!grepl("_fit", names(r$fits))]
   
    perf <- get_performance_table(gs, d)
    
    # remember for easy plotting
    perf <- mutate(
      perf,
      rdegree = as.character(r$rdegree),
      snp = r$snp,
      iteration = iteration,
      density_true =
        edge_density(igraph.from.graphNEL(r$graph.observed))
    )
    perf
  }) %>% bind_rows() # lapply over simulations ---------------------------------
  
  return(tab)
}

#'------------------------------------------------------------------------------
#' Use the bdgraph internal method to get spec/sens, f1 and MCC. Uses the
#' original simulation object containing the ground truth graph
#' 
#' @param graph_list Named list of graph objects
#' @param data_sim BDgraph data simulation object from which graphs have been 
#' inferred
#' 
get_performance_table <- function(graph_list, data_sim) {
  
  require(dplyr)
  require(igraph)
  
  perf <- lapply(names(graph_list), function(g) {

    current_graph <- graph_list[[g]]
    if(is.null(current_graph)) {
    
      return(NULL)

    }

    # we need the adjacency matrix for comparison
    perf <- t(BDgraph::compare(data_sim, as(current_graph, "matrix")))
    comparisons <- c("True", g)
    perf <- as.data.frame(perf)
    rownames(perf) <- paste(g, comparisons, sep = "_")
    perf <- perf[!grepl("True", rownames(perf)), ]
    
    # annotate density for comparison
    ig <- igraph::igraph.from.graphNEL(current_graph)
    dens <- edge_density(ig)
    perf$density_model <- dens
    perf %>% mutate(comparison = g)
    
  }) %>% bind_rows()
  
  return(perf)
}
# ------------------------------------------------------------------------------
#' Generate the validation table for all simul results contained in the provided
#' data object.
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' 
# ------------------------------------------------------------------------------
get_validation_table_prior_completeness <- function(result, iteration) {
  
  require(dplyr)
  require(igraph)
  
  results <- result$fits
  d <- result$data.sim
  
  # check each individual simulated result
  tab <- lapply(names(results), function(n) {
    
    r <- results[[n]]
    
    # Retain only the graph objects
    gs <- r[!grepl("_fit", names(r))]
    
    perf <- get_performance_table(gs, d)
    
    # remember more information for easier plotting
    perf <- mutate(
      perf,
      fraction_to_keep = n,
      snp = result$snp,
      iteration = iteration,
      comparison = names(gs),
      density_true =
        edge_density(igraph.from.graphNEL(result$graph.observed))
    )
    perf
  }) %>% bind_rows() # lapply over different fractions
  
  return(tab)
}

# ------------------------------------------------------------------------------
#' Run the GGM inference for all models on all simulated data (i.e. for
#' different noise degrees)
#' 
#' @param simulated_data The main simulation object containing all simulated
#' data and graphs for differing prior noise degrees
#' @param priors The prior matrix. Will not be used for the rbinom prior which
#' is provided in the respected simulation object
#' @param ranges The original ranges collection for the related locus
#' @param fcpg_context The TF-cpg annotation context (for graph annotation)
#' @param ppi_db The underlying PPI network (for graph annotation)
#' @param threads The number of threads which can be used, Default: 1
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#' 
# ------------------------------------------------------------------------------
run_ggm <- function(simulated_data, priors, ranges, 
                    fcpg_context, ppi_db, 
                    subset = c("all", "minimal", seq(50,700,by=50)),
                    threads = 1) {
  
  # check subset. by default we use all data. in case a numeric subset is
  # specifified, we only run the inference for the non-noisy prior model.
  subset <- match.arg(subset)

  print("Subset is:")
  print(subset)
  
  # iterate over all simulations (different noise levels)
  # in case subset != all, we only process the non-noisy prior
  sims <- names(simulated_data)
  if(subset != "all") {
    
    # special case: subset + iterate over all noise levels
    if(subset == "minimal") {
      subset <- as.numeric(snakemake@params$minimal_subset_size)     
    } else {
      sims <- sims[grepl("rd0$|rd0.8$", sims)]
      subset <- as.numeric(subset)
    }
  }
  
  result <- lapply(sims, function(n) {
    # --------------------------------------------------------------------------
    # Get data
    
    sim <- simulated_data[[n]]
    # sentinel name and simulated data
    s <- sim$snp
    d <- sim$data.sim$data

    # in case of rbinom simulation, we get adjusted priors
    if (grepl("_rbinom", n)) {
      priors <- sim$priors
    }
    
    # check whether we need a subset
    if(is.numeric(subset)) {
        d <- d[sample(1:nrow(d), size=subset, replace = F),]
    }
    
    # 
    # --------------------------------------------------------------------------
    print("Infer regulatory networks.")
    
    result <-
      infer_all_graphs(d, priors, ranges, fcpg_context, ppi_db, threads)
    sim$fits <- result
    
    sim
  })
  names(result) <- sims
  
  return(result)
}


# ------------------------------------------------------------------------------
#' Run the GGM inference for all models on simulated data for different 'levels'
#' of prior completeness (based on quantiles)
#' 
#' @param simulated_data The main simulation object containing all simulated
#' data and graphs for differing prior noise degrees
#' @param priors The prior matrix. Will not be used for the rbinom prior which
#' is provided in the respected simulation object
#' @param ranges The original ranges collection for the related locus
#' @param fcpg_context The TF-cpg annotation context (for graph annotation)
#' @param ppi_db The underlying PPI network (for graph annotation)
#' @param threads The number of threads which can be used, Default: 1
#' 
#' @author Johann Hawe <johann.hawe@tum.de>
#' 
# ------------------------------------------------------------------------------
run_ggm_prior_completeness <- function(simulated_data, priors, ranges, 
                                       fcpg_context, ppi_db, threads = 1) {
 
  # use only the 'zero' noise level data
  sims <- names(simulated_data)
  sim <- simulated_data[[sims[grepl("rd0$", sims)]]]
  
  fractions <- seq(0.1,0.9, by=0.1)
  
  results <- lapply(fractions, function(fraction_to_keep) {
    
    print(paste0("Current fraction: ", fraction_to_keep))
    
    # sentinel name and simulated data -----------------------------------------
    s <- sim$snp
    d <- sim$data.sim$data
    
    # Drop/adjust priors -------------------------------------------------------
    pseudo_prior <- min(priors)
    total_prior_edges <- sum(priors[upper.tri(priors)] > pseudo_prior)
    priors <- drop_prior_edges(priors, fraction_to_keep, pseudo_prior)
    print("Fraction of priors remaining:")
    print(sum(priors[upper.tri(priors)] > pseudo_prior) / total_prior_edges)
     
    # --------------------------------------------------------------------------
    print("Infer regulatory networks.")
    
    result <-
      infer_all_graphs(d, priors, ranges, fcpg_context, ppi_db, threads)

    result
  })
  names(results) <- paste0("fraction", fractions)
  sim$fits <- results
  
  return(sim)
}

#' -----------------------------------------------------------------------------
#' Drops a certain fraction of prior edges from a prior matrix. Edges are set
#' to the provided pseudo_prior
#' 
#' @param priors The prior matrix
#' @param fraction_to_keep The fraction of prior edges to keep (e.g. 0.2)
#' @param pseudo_prior Value of the pseudo prior to set for 'dropped' edges (e.g. 1e-7)
#' -----------------------------------------------------------------------------
drop_prior_edges <- function(priors,
                             fraction_to_keep, 
                             pseudo_prior) {
  
  # temp save priors for melting
  p <- priors
  p[upper.tri(p, TRUE)] <- NA
  p <- melt(p, na.rm = T)
  
  idxs <- which(p$value > pseudo_prior)
  total_priors <- length(idxs)
  total_to_keep <- ceiling(total_priors * fraction_to_keep)
  total_to_drop <- total_priors - total_to_keep
  
  set.seed(42)
  
  idxs_to_drop <- idxs[!sample(c(rep(TRUE, total_to_keep),
                                 rep(FALSE, total_to_drop)))]
  
  # set pseudo prior for all entries
  for (i in idxs_to_drop) {
    n1 <- as.character(p[i, "Var1"])
    n2 <- as.character(p[i, "Var2"])
    priors[n1, n2] <- priors[n2, n1] <- pseudo_prior
  }
  
  return(priors)
}
