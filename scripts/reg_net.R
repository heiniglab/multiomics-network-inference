#' -----------------------------------------------------------------------------
#' The main script handling the inference of regulatory network from genomic
#' data.
#' The main method "reg_net()" allows to infer networks using several distinct
#' approaches, for some also providing different default parameters.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Dec 13 17:49:42 2018
#' -----------------------------------------------------------------------------

# when sourcing, check which packages are installed for quicker diagnostics
require(BDgraph)
require(glasso)
require(GENIE3)
require(iRafNet)
require(GeneNet)
require(doParallel)
require(parallel)
require(cvTools)
require(igraph)
require(graph)
require(reshape2)

# define the available models for network inference
reg_net.models <- function() {
  return(c("genenet", "bdgraph", "irafnet", "glasso", "genie3" , "custom"))
}

# ------------------------------------------------------------------------------
#' Main method for inferring regulatory networks from different models
#'
#' @param data The data matrix (n x p) from which to infer the network. Must
#' have column names (i.e. node names)
#' @param priors The matrix of priors (p x p) which to use. Needs to be set to
#' NULL explicitely to disregard priros. In case of bdgraph model, these are
#' used to define a start graph for the algorithm.
#' @param model The model to be used. One of `reg_net.models()`
#' @param threads optional. Number of threads to be used (if applicable)
#' @param use_gstart optional. Flag whether to use prior based start
#' configuration in case of BDgraph Default: TRUE
#' @param iter optional. Number of iterations to be performed for BDgraph
#' @param burnin optional. Number of burnin iterations for BDgraph
#' @param ntrees optional. Number of trees to build per variable in iRafNet
#' @param mtry optional. Number of variables to select from forests in iRafNet.
#' Default: round(sqrt(ncol(data)-1))
#' @param npermut optional. Number of permutations to perform for
#' iRafNet background
#'
#' @return Returns a list containing bot the final model fit as well as the
#' graphNEL object extracted from that model fit.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
# ------------------------------------------------------------------------------
reg_net <- function(data,
                    priors,
                    model,
                    threads = 1,
                    use_gstart = T,
                    gstart = NULL,
                    iter = 10000,
                    burnin = 5000,
                    ntrees = 1000,
                    mtry = round(sqrt(ncol(data) - 1)),
                    npermut = 5,
                    verbose = FALSE) {
  require(doParallel)
  
  # vector of all nodes
  nodes <- colnames(data)
  
  # get available models
  ms <- reg_net.models()
  if (!(model %in% ms)) {
    stop(paste0("Model not supported: ", model))
  } else {
    print(paste0("Running model '", model, "'"))
  }
  
  if (any(is.na(data))) {
    warning("Handling NAs present in data automatically!")
  }
  
  # for genenet and iRafNet, remove NAs
  data_no_nas <- data[, apply(data, 2, function(x)
    ! anyNA(x))]
  
  if (!is.null(priors)) {
    # if prior dimensions do not match, try to subset priors
    if (ncol(data) < ncol(priors)) {
      warning("Subsetting priors automatically. You should check your input.")
      priors <- priors[colnames(data), colnames(data)]
      
      # sanity check: do colnames match?
      if (!all(colnames(priors) == colnames(data))) {
        stop("Column names did not match after autmatic prior subsetting.")
      }
    }
    # we possibly lost some data when filtering for NAs, so we adjust priors too
    priors_no_nas <-
      priors[rownames(priors) %in% colnames(data_no_nas),
             colnames(priors) %in% colnames(data_no_nas)]
  } else {
    priors_no_nas <- NULL
  }
  # check which model to build
  if ("genenet" %in% model) {
    require(GeneNet)
    
    # get partial correlation estimates
    pcors <- ggm.estimate.pcor(data.matrix(data_no_nas), verbose = FALSE)
    
    fit <- network.test.edges(pcors, plot = FALSE, verbose = FALSE)
    fit$node1 <- colnames(data_no_nas)[fit$node1]
    fit$node2 <- colnames(data_no_nas)[fit$node2]
    class(fit) <- c(class(fit), "genenet")
    
  } else if ("bdgraph" %in% model) {
    require(BDgraph)
    # check some conditions before applying the model
    # no priors
    if (is.null(priors)) {
      if(verbose) {
        print("Setting uniform prior.")
      }
      
      # set default uniform prior
      # 0.5 reflects uniform according to email from the package author from 4.8.17
      bdpriors <- 0.5
      
      # build without priors, nevertheless, we could have a custom start graph
      if (use_gstart) {
        if (is.null(gstart)) {
          stop("Start graph must be provided if no priors are given.")
        }
      } else {
        gstart <- "empty"
      }
    } else {
      # we have priors
      bdpriors <- priors
      
      # check whether to use the prior based start graph
      if (use_gstart) {
        gstart <- get_gstart_from_priors(bdpriors)
      } else {
        gstart <- "empty"
      }
    }
    # perform the model fitting
    fit <- try({
      bdgraph(
        data,
        method = "gcgm",
        g.prior = bdpriors,
        iter = iter,
        burnin = burnin,
        g.start = gstart,
        save = T,
        cores = threads
      )
    })
    
  } else if ("irafnet" %in% model) {
    require(iRafNet)
    
    
    cl <- makeCluster(threads)
    registerDoParallel(cl)
    
    irn_out <-
      iRafNet(data_no_nas,
              priors_no_nas,
              ntrees,
              mtry,
              colnames(data_no_nas),
              threads = threads)
    
    irn_perm_out <- Run_permutation(
      data_no_nas,
      priors_no_nas,
      ntrees,
      mtry,
      colnames(data_no_nas),
      npermut,
      threads = threads
    )
    
    fit <- list(
      irn_out = irn_out,
      irn_perm_out = irn_perm_out,
      nodes = colnames(data_no_nas)
    )
    
    class(fit) <- c(class(fit), "irafnet")
  
    stopCluster(cl)
    
  } else if ("glasso" %in% model) {
    require(glasso)
    if (!is.null(priors)) {
      gl_out <- glasso_cv(data, priors, nodes, 
                          threads = threads, verbose = verbose)
    } else {
      gl_out <- glasso_cv(data, NULL, nodes, 
                          threads = threads, verbose = verbose)
    }
    class(gl_out) <- c(class(gl_out), "glasso")
    fit <- gl_out
  } else if ("genie3" %in% model) {
    require(GENIE3)
    fit <- genie3(data_no_nas, 
                  threads = threads, verbose = verbose)
    class(fit) <- c(class(fit), "genie3")
  } else if ("custom" %in% model) {
    stop("Sorry, custom model is not yet implemented.")
  }
  
  # now get the graph object
  g <- graph_from_fit(fit, nodes, annotate = F, verbose = verbose)
  
  return(list(graph = g, fit = fit))
}

# ------------------------------------------------------------------------------
#' Creates a graphNEL object from a given model result
#' Evaluates fit to powerlaw distribution to determine significance cutoffs
#' for some models (currently: genenet, genie3 and irafnet)
#'
#' @param ggm.fit The ggm fit
#' @param nodes All nodes of the graph. Will be used to create the
#' graph object
#' @param ranges The ranges of the entities used for the graph fit
#' @param string_db
#' @param annotate Flag whether to annotate the graph entities with nodeData.
#' Default: T
#'
#' @value graph-nel object created from the bdgraph result
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
graph_from_fit <- function(ggm.fit,
                           nodes,
                           ranges = NULL,
                           ppi_db = NULL,
                           annotate = T, 
                           verbose = FALSE) {
  if (annotate & is.null(ranges)) {
    stop("Ranges must not be null for annotating graphs!")
  }
  
  library(BDgraph)
  library(graph)
  library(igraph)
  require(reshape2)
  
  # get the graph instance from the ggm fits
  if (inherits(ggm.fit, "bdgraph")) {
    best_cutoff <-
      get_best_graph_cutoff(
        seq(0.5, 0.99, by = 0.01),
        get_bdgraph_graph,
        threads = 1,
        rsquare.cut = 0.8,
        ggm.fit,
        nodes
      )
    
    if(verbose) {
      print(paste0("BDgraph best cutoff is: ", best_cutoff))
    }
    
    g <- get_bdgraph_graph(best_cutoff, ggm.fit, nodes)
    
  } else if (inherits(ggm.fit, "irafnet")) {
    best_cutoff <-
      get_best_graph_cutoff(
        seq(0.01, 0.2, by = 0.01),
        get_irafnet_graph,
        threads = 1,
        rsquare.cut = 0.8,
        ggm.fit$irn_out,
        ggm.fit$irn_perm_out,
        nodes
      )
    
    g <-
      get_irafnet_graph(best_cutoff, ggm.fit$irn_out, ggm.fit$irn_perm_out,
                        nodes)
    
  } else if (inherits(ggm.fit, "genenet")) {
    best_cutoff <-
      get_best_graph_cutoff(
        seq(0.8, 0.99, by = 0.01),
        get_genenet_graph,
        threads = 1,
        rsquare.cut = 0.8,
        ggm.fit,
        nodes
      )
    g <- get_genenet_graph(best_cutoff, ggm.fit, nodes)
    
  } else if (inherits(ggm.fit, "glasso")) {
    pm <- ggm.fit$wi
    cn <- colnames(pm)
    g <- graphNEL(nodes, edgemode = "undirected")
    
    # create edge matrix
    v <- sm2vec(pm)
    vidx <- sm.index(pm)
    em <- cbind.data.frame(
      pcor = v,
      node1 = cn[vidx[, 1]],
      node2 = cn[vidx[, 2]],
      stringsAsFactors = F
    )
    
    # define edges with pcor != 0
    em <- subset(em, pcor != 0)
    if (nrow(em) > 0) {
      g <- addEdge(em$node1, em$node2, g)
    }
  } else if (inherits(ggm.fit, "genie3")) {
    g <- get_genie3_graph(ggm.fit$best_weight,
                          nodes,
                          ggm.fit$linklist)
  }
  
  if (annotate) {
    # set node and edge attributes
    g <- annotate.graph(g, ranges, ppi_db)
  }
  return(g)
}

# ------------------------------------------------------------------------------
# Gets a graph object based on a BDgraph graph fit and a specific posterior
# probability cutoff
# ------------------------------------------------------------------------------
get_bdgraph_graph <- function(posterior_probability,
                              model_fit,
                              nodes) {
  # only gives upper.tri -> invert
  g_adj <- t(BDgraph::select(model_fit, cut = posterior_probability))
  cn <- colnames(g_adj)
  g <- graphNEL(nodes,
                edgemode = "undirected")
  
  # create edge matrix
  v <- sm2vec(g_adj)
  vidx <- sm.index(g_adj)
  em <- cbind.data.frame(
    is_edge = v,
    node1 = cn[vidx[, 1]],
    node2 = cn[vidx[, 2]],
    stringsAsFactors = F
  )
  # get edges
  em <- subset(em, is_edge == 1)
  if (nrow(em) > 0) {
    g <- addEdge(em$node1, em$node2, g)
  }
  return(g)
}

# ------------------------------------------------------------------------------
#' Gets a graph object based on a iRafNet graph fit and a specific FDR cutoff
# ------------------------------------------------------------------------------
get_irafnet_graph <- function(fdr_cutoff,
                              model_output,
                              permutation_output,
                              nodes) {
  fit <- iRafNet_network(model_output,
                         permutation_output,
                         TH = fdr_cutoff)
  
  g <- graphNEL(nodes, edgemode = "undirected")
  g <- addEdge(fit$gene1, fit$gene2, g)
  
  return(g)
}

# ------------------------------------------------------------------------------
#' Gets a graph object from a genenet model fit and given a specific link
#' probability cutoff (1-fdr)
# ------------------------------------------------------------------------------
get_genenet_graph <- function(prob_cutoff, model_fit, nodes) {
  net <- extract.network(model_fit,
                         cutoff.ggm = prob_cutoff,
                         verbose = FALSE)
  g <- graphNEL(nodes, edgemode = "undirected")
  
  if (nrow(net) > 0) {
    g <- addEdge(net$node1, net$node2, g)
  }
  return(g)
}

# ------------------------------------------------------------------------------
#' Get a graph score summary to estimate how well the inferred graph structure
#' fits our assumptions of the underlying regulatory mechanisms.
#'
#' @param g The graph for which to get the score (igraph)
#' @param sentinel The sentinel to be found in the graph
#' @param ranges The ranges collection for the respective sentinel/locus
#' @param density optional. Used to adjust final score (the higher the density,
#' the lower the score)
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
get_graph_score <- function(g, sentinel, ranges, density = NULL) {
  require(igraph)
  
  seed <- ranges$seed
  
  # get graph score
  score <- 0
  
  # keep only cluster with sentinel
  cl <- clusters(g)
  if (!sentinel %in% names(cl$membership))
    return(0)
  
  cn <- cl$membership[sentinel]
  subnodes <- names(cl$membership[cl$membership == cn])
  if (length(subnodes) < 2)
    return(0)
  
  g <- induced_subgraph(g, subnodes)
  
  # check whether we have any SNP gene.
  # if not, the score will be 0
  # otherwise proceed with check
  adj_sent <- igraph::neighbors(g, sentinel)$name
  
  # get number of trans genes in our cluster
  v <- V(g)$name
  if (seed == "meqtl") {
    tg_in_v <- intersect(names(ranges$cpgs), v)
  } else {
    tg_in_v <- intersect(ranges$trans_genes$SYMBOL, v)
  }
  
  if (seed %in% c("meqtl", "eqtlgen", "eqtl")) {
    # weight for cis-genes
    X <- 0.5
    cis_sent <- intersect(ranges$snp_genes$SYMBOL, adj_sent)
    if (length(cis_sent) < 1) {
      return(score)
    } else {
      score <- score + X
    }
    
    # define proportion to use for scoring
    cis_in_v <- intersect(ranges$snp_genes$SYMBOL, v)
    cis_not_sent <- setdiff(cis_in_v, cis_sent)
    Y_cis <- X / length(cis_in_v)
    Y_trans <- X / length(tg_in_v)
    
    # check for cis genes which are not connected to the snp/any other
    # snp genes
    for (cis in cis_not_sent) {
      # check adjacency to any of the cis genes adj to sentinel
      if (!any(cis_sent %in% neighbors(g, cis))) {
        score <- score - Y_cis
      }
    }
    
    # check whether trans genes are directly connected to sentinel
    # score gets a penalty for each such instance
    tg_sent <- intersect(tg_in_v, adj_sent)
    for (i in tg_sent) {
      score <- score - Y_trans
    }
    
    # now we investigate whether we can reach the trans genes
    # starting from the snp_genes. for this we first
    # remove the sentinel from the network
    g <- delete.vertices(g, sentinel)
    for (cis in cis_sent) {
      for (trans in tg_in_v) {
        # we always remove all other trans genes than the current
        # one to avoid walking over other trans genes to reach it
        sg <- delete.vertices(g, setdiff(tg_in_v, trans))
        
        # get path to trans_gene
        paths <-
          suppressWarnings(get.shortest.paths(sg, cis, trans))$vpath[[1]]
        if (length(paths) > 1) {
          score <- score + Y_trans
        }
      }
    }
  } else {
    # TODO: we could check whether sentinel is TF, then get fraction of
    # directly reachable trans genes
    
    # fraction of reachable TGs via at least one TF or spath gene
    Y_trans <- 1 / length(tg_in_v)
    
    for (trans in tg_in_v) {
      # remove all other trans-genes to get only independently reachable TGs
      sg <- delete.vertices(g, setdiff(tg_in_v, trans))
      
      paths <-
        suppressWarnings(get.shortest.paths(sg, sentinel, trans))$vpath[[1]]
      
      # we need at least one TF or Spath gene on the path
      pgenes <- unique(c(ranges$tfs$SYMBOL, ranges$spath$SYMBOL))
      
      if (length(paths) > 1 && any(pgenes %in% paths$name)) {
        score <- score + Y_trans
      }
    }
  }
  # if graph density is available, we use it to adjust the score
  if (!is.null(density)) {
    # TODO might want to adjust density weighting, e.g. score*(-log10(density))
    score <- score * (-log10(density))
  }
  return(score)
}


#' -----------------------------------------------------------------------------
#' Calculates loglikelihood according to
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3413177/
#'
#' @param theta The estimated precision matrix (e.g. on training data)
#' @param sigma The sample covariance matrix (e.g. on test data)
#' @param n The number of samples
#' @param p The number of parameters
#'
#' -----------------------------------------------------------------------------
loglik <- function(theta, sigma, n, p) {
  c1 <- 0.5 * n
  c2 <- p * log(2 * pi)
  logdet <- determinant(theta, logarithm = T)$modulus
  ll <- c1 * (logdet - sum(sigma * theta) - c2)
  
  return(ll)
}

#' -----------------------------------------------------------------------------
#' Calculate the BIC for a specific model parameters.
#'
#' @param ll The log likelihood
#' @param n The number of samples
#' @param k The number of estimated parameters
#'
#' -----------------------------------------------------------------------------
bic <- function(ll, n, k) {
  return(log(n) * k - 2 * ll)
}

#' -----------------------------------------------------------------------------
#' Perform glasso cross validation for lambda selection
#'
#' Allows the use of a prior matrix (argument priors) which will be supplied
#' to the model fitting if provided. Will be scaled by the Rho.
#'
#' NOTE: code was adapted from
#' https://github.com/czarrar/causality_primer/blob/master/journal/2014-07-12_cross_validation_graphical_lasso.Rmd
#'
#' @param data the data matrix for which to fit the model (n x p)
#' @param priors The prior matrix (p x p). Default: NULL
#' @param k The number of CV folds to use
#' @param rholist The list of rhos to check
#' @param threads the number of threads which can be used
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
glasso_cv <- function(data,
                      priors = NULL,
                      nodes,
                      k = 5,
                      rholist = seq(0.01, 1, by = 0.005),
                      threads = 1,
                      verbose = FALSE) {
  require(glasso)
  require(cvTools)
  require(parallel)
  
  if(verbose) {
    print("Starting glasso CV.")
  }
  
  n <- nrow(data)
  folds <- cvFolds(n, k, type = "random")
  
  # define penalties to be used in glasso
  if (!is.null(priors)) {
    penalties <- (1 - priors)
  }
  
  # get the ll progression for all CV folds
  loglikes <- lapply(1:k, function(ki) {
    S_train <- cov(data[folds$which != ki,], use = "na.or.complete")
    S_test  <-
      cov(data[folds$which == ki,], use = "na.or.complete")
    
    # get the number of variables
    p <- ncol(S_test)
    n <- nrow(S_test)
    
    # check all rhos
    # NOTE: we can't use glassopath since we want to include the prior
    # information
    unlist(mclapply(rholist, function(rho) {
      if (!is.null(priors)) {
        gl <- glasso(S_train,
                     rho = rho * penalties,
                     penalize.diagonal = F)
      } else {
        gl <- glasso(S_train,
                     rho = rho,
                     penalize.diagonal = F)
      }
      rownames(gl$wi) <- colnames(gl$wi) <- colnames(data)
      class(gl) <- c(class(gl), "glasso")
      
      # get log likelihood, number of edges and BIC
      ll <- loglik(gl$wi, S_test, n, p)
      
      g <- graph_from_fit(gl, nodes, annotate = F)
      ne <- numEdges(g)
      
      # eBIC?
      bic <- bic(ll = ll, n = n, k = ne)
      
      c(ll, bic, ne)
      
    }, mc.cores = threads))
  })
  loglikes <- do.call(rbind, loglikes)
  colnames(loglikes) <-
    paste0(c("rho_rho=", "bic_rho=", "edgeNum_rho="),
           rep(rholist, each = 3))
  rownames(loglikes) <- paste0("k=", 1:k)
  
  # get the 'best' rho according to the BIC
  bic_sub <- colMeans(loglikes[, grepl("bic", colnames(loglikes))])
  ind <- which.min(bic_sub)
  
  # retrieve rho from the name
  rho <- as.numeric(strsplit(names(ind), "=")[[1]][2])
  
  # retrain the full model
  S <- cov(data, use = "na.or.complete")
  if (!is.null(priors)) {
    model <- glasso(S, rho * penalties, penalize.diagonal = F)
  } else {
    model <- glasso(S, rho, penalize.diagonal = F)
  }
  # set names on precision matrix
  rownames(model$wi) <- colnames(model$wi) <- colnames(data)
  
  # remember rho and loglikelihoods
  model$rho_best <- rho
  model$rho_progression <- loglikes
  
  return(model)
}

#' -----------------------------------------------------------------------------
#' Gets a graphnel object from a fitted GENIE3 model for a specific link weight
#' threshold
#'
#' @param threshold The treshold to be applied in GENIE3's \code{getLinkList}
#' @param nodes All nodes used in the calculation of the model
#' @param linklist The ranked list of links as extracted from GENIE3
#'
#' -----------------------------------------------------------------------------
get_genie3_graph <- function(threshold, nodes, linklist) {
  require(graph)
  
  edges <- subset(linklist, weight >= threshold)
  
  if (nrow(edges) > 0) {
    g <- graph::graphNEL(nodes,
                         edgemode = "undirected")
    g <- graph::addEdge(as.character(edges$regulatoryGene),
                        as.character(edges$targetGene),
                        g)
  } else {
    g <- graph::graphNEL(nodes, edgemode = "undirected")
  }
  return(g)
}

#' -----------------------------------------------------------------------------
#' Trains an optimized GENIE3 model
#'
#' Uses GENIE3 to estimate a regulatory network from the supplied data matrix.
#' The threshold for the regulatory links is optimized by fitting the graph
#' for a specific threshold to a powerlaw and then selecting the threshold for
#' which the likelihood of the fitted model is maximal
#'
#' @param data The data matrix (n x p)
#' @param threads The number of threads to be used
#' @param verbose Whether to be verbose about the process or not
#'
#' -----------------------------------------------------------------------------
genie3 <-
  function(data,
           threads = 1,
           verbose = FALSE) {
    require(GENIE3)
    require(igraph)
    require(parallel)
    
    # we expect a n x p matrix, genie3 needs a p x n matrix
    if (verbose) {
      start <- Sys.time()
    }
    model <- GENIE3(t(data), nCores = threads)
    linklist <- getLinkList(model, threshold = 0)
    if (verbose) {
      end <- Sys.time()
      print(paste0("GENIE3 benchmark: ", end - start))
    }
    all_link_weights <- sort(unique(linklist$weight))
    
    lw <- length(all_link_weights)
    # get a subset to check if needed
    if (lw > 500) {
      all_link_weights <- quantile(all_link_weights,
                                   seq(0, 1, by = 0.005))
    }
    
    if (verbose) {
      print(paste0("Checking ", length(all_link_weights), " weights."))
    }
    best_weight <-
      get_best_graph_cutoff(
        all_link_weights,
        get_genie3_graph,
        threads = 1,
        rsquare.cut = 0.8,
        colnames(data),
        linklist
      )
    
    return(
      list(
        model = model,
        nodes = colnames(data),
        linklist = linklist,
        best_weight = best_weight
      )
    )
  }

# ------------------------------------------------------------------------------
#' Gets the best cutoff from a list of given cutoffs with respect to the
#' goodness of fit to a powerlaw distribution of graphs
#'
#' @param cutoff_list The list of cutoff values to be screened (as a vector)
#' @param model_type The type of the model to get the best cutoff for. Currently
#' either 'genie3', 'irafnet' or 'genenet'
#' @param graph_extraction_callback A method which can be used to obtain a new
#' graph object for the specified model. First param needs to be the currently selected cutoff/weight.
#' @parma threads Number of threads to be used for screening
#' @param ... Extra parameteres being passed to the graph extraction callback
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
get_best_graph_cutoff <- function(cutoff_list,
                                  graph_extraction_callback,
                                  threads = 1,
                                  rsquare.cut = 0.8,
                                  ...) {
  powerlaw_fits <- mclapply(cutoff_list, function(cutoff) {
    g <- graph_extraction_callback(cutoff, ...)
    fit <- get_powerlaw_fit(g)
    if (!is.null(fit)) {
      fit$weight <- cutoff
    }
    fit
  }, mc.cores = threads)
  
  powerlaw_fits <- do.call(rbind.data.frame, powerlaw_fits)
  colnames(powerlaw_fits) <-
    c("ll", "ks_p", "alpha", "beta", "r2", "mcon", "weight")
  rownames(powerlaw_fits) <- paste0("weight=", powerlaw_fits$weight)
  
  # get best weight (r2>=cutoff, then highest mean connectivity)
  fits_r2 <- subset(powerlaw_fits, r2 >= rsquare.cut)
  if (nrow(fits_r2) < 1) {
    # didn't find a cutoff, so we take the maximum r2
    warning(paste0("No R^2 above ", rsquare.cut))
    fits_r2 <- subset(powerlaw_fits, r2 == max(r2))
  }
  fits_mcon <- subset(fits_r2, mcon == max(mcon))
  fits_beta <- fits_mcon[which.min(-1 - fits_mcon$beta), ]
  
  return(fits_beta$weight)
}

#' -----------------------------------------------------------------------------
#' we can use igraph packages 'fit_power_law' to fit the degree distribution
#' to a power law distr for each weight cutoff
#'
#' @param graph The graph object for which to get the fit
#' @param nbreaks The number of breaks to use when discretizing the degree
#' distribution / getting its frequencies. Default: 20
#'
#' @author Johann Hawe <johann.hawe@tum.de>
#' -----------------------------------------------------------------------------
get_powerlaw_fit <- function(graph, nbreaks = 20) {
  ds <- graph::degree(graph)
  
  if (var(ds) == 0)
    return(NULL)
  
  ds1 <- ds + 1
  # fit power law
  fit <- igraph::fit_power_law(ds1)
  
  # get the correlation between log10(p(k)) and log10(k)
  # taken from: https://github.com/cran/WGCNA/blob/master/R/Functions.R
  # function: scaleFreeFitIndex
  
  # devide degrees into bins and get the respective frequencies
  discretized = cut(ds, nbreaks)
  dk = as.vector(tapply(ds, discretized, mean))
  pdk = as.vector(tapply(ds, discretized, length) / length(ds))
  pdk = ifelse(is.na(pdk), 0, pdk)
  
  breaks1 = seq(from = min(ds),
                to = max(ds),
                length = nbreaks + 1)
  hist1 = hist(ds,
               breaks = breaks1,
               plot = FALSE,
               right = TRUE)
  dk2 = hist1$mids
  dk = ifelse(is.na(dk), dk2, dk)
  dk = ifelse(dk == 0, dk2, dk)
  
  # get logs of degree distr and frequency, then calculate their linear fits
  log_dk = as.vector(log10(dk))
  log_pdk = as.numeric(log10(pdk + 1e-09))
  lm1 = lm(log_pdk ~ log_dk)
  r2 <- summary(lm1)$r.squared
  beta <- lm1$coefficients["log_dk"]
  
  # we need a negative slope for a 'plausible' network
  if (beta > 0)
    r2 <- r2 * (-1)
  
  # get the mean node connectivity
  mcon <- (sum(ds - 1)) / 2
  
  # return ll of fitted distr
  list(
    ll = fit$logLik,
    ks_p = fit$KS.p,
    alpha = fit$alpha,
    beta = beta,
    r2 = r2,
    mcon = mcon
  )
}

#' -----------------------------------------------------------------------------
#' Infers graphs for all available models and annotates them.
#'
#' @param data data matrix (n x p) to be used for fitting the models
#' @param priors the p x p prior matrix
#' @param ranges the ranges collection underlying the data
#' @param fcontext the tfbs context file
#' @param ppi_db the ppi db (graph) underlying the data
#' @param threads number of threads to be used. Default: 1
#'
#' -----------------------------------------------------------------------------
infer_all_graphs <-
  function(data,
           priors,
           ranges,
           fcontext,
           ppi_db,
           threads = 1) {
    
    # we set the OMP/BLAS number of threads to 1
    # this avoids issues we had in the glasso CV with multi-threading on cluster
    # also necessary for BDgraph
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    print("Fitting model using GeneNet.")
    genenet <- reg_net(data, NULL, "genenet", threads = threads)
    
    print("Fitting model using glasso.")
    glasso <- reg_net(data, priors, "glasso", threads = threads)
    
    print("Fitting model using glasso, no priors.")
    glasso_no_priors <-
      reg_net(data, NULL, "glasso", threads = threads)
    
    print("Fitting model using genie3.")
    genie3 <- reg_net(data, NULL, "genie3", threads = threads)
    
    print("Fitting bdgraph using priors.")
    bdgraph <- reg_net(data, priors, "bdgraph", threads = threads)
    
    print("Fitting bdgraph without priors.")
    bdgraph_no_priors <- reg_net(data,
                                 NULL,
                                 "bdgraph",
                                 use_gstart = F,
                                 threads = threads)
    
    print("Fitting bdgraph without priors, full start graph.")
    bdgraph_no_priors_full <- reg_net(
      data,
      NULL,
      "bdgraph",
      use_gstart = T,
      gstart = "full",
      threads = threads
    )
    
    print("Fitting model using iRafNet.")
    irafnet <- reg_net(data, priors, "irafnet", threads = threads)
    
    # --------------------------------------------------------------------------
    print("Add custom annotations for the graphs.")
    # --------------------------------------------------------------------------
    bdgraph$graph <- annotate.graph(bdgraph$graph,
                                    ranges, ppi_db, fcontext)
    
    bdgraph_no_priors$graph <-
      annotate.graph(bdgraph_no_priors$graph,
                     ranges, ppi_db, fcontext)
    
    bdgraph_no_priors_full$graph <-
      annotate.graph(bdgraph_no_priors_full$graph,
                     ranges, ppi_db, fcontext)
    
    irafnet$graph <- annotate.graph(irafnet$graph,
                                    ranges, ppi_db, fcontext)
    
    genenet$graph <- annotate.graph(genenet$graph,
                                    ranges, ppi_db, fcontext)
    
    genie3$graph <- annotate.graph(genie3$graph,
                                   ranges, ppi_db, fcontext)
    
    glasso$graph <- annotate.graph(glasso$graph,
                                   ranges, ppi_db, fcontext)
    
    glasso_no_priors$graph <- annotate.graph(glasso_no_priors$graph,
                                             ranges, ppi_db, fcontext)
    
    # --------------------------------------------------------------------------
    print("Create result list.")
    # --------------------------------------------------------------------------
    result <- list(
      # bdgraph
      bdgraph_fit = bdgraph$fit,
      bdgraph = bdgraph$graph,
      # bdgraph no priors
      bdgraph_no_priors_fit = bdgraph_no_priors$fit,
      bdgraph_no_priors = bdgraph_no_priors$graph,
      # bdgraph no priors, full start
      bdgraph_no_priors_full_fit = bdgraph_no_priors_full$fit,
      bdgraph_no_priors_full = bdgraph_no_priors_full$graph,
      # irafnet
      irafnet_fit = irafnet$fit,
      irafnet = irafnet$graph,
      # genenet
      genenet_fit = genenet$fit,
      genenet = genenet$graph,
      # glasso
      glasso_fit = glasso$fit,
      glasso = glasso$graph,
      # glasso no priors
      glasso_no_priors_fit = glasso_no_priors$fit,
      glasso_no_priors = glasso_no_priors$graph,
      # genie3
      genie3_fit = genie3$fit,
      genie3 = genie3$graph
    )
    
    return(result)
  }


#' -----------------------------------------------------------------------------
#' Infers graphs for all available prior based models and annotates them.
#'
#' @param data data matrix (n x p) to be used for fitting the models
#' @param priors the p x p prior matrix
#' @param ranges the ranges collection underlying the data
#' @param fcontext the tfbs context file
#' @param ppi_db the ppi db (graph) underlying the data
#' @param threads number of threads to be used. Default: 1
#'
#' -----------------------------------------------------------------------------
infer_all_graphs_priors <-
  function(data,
           priors,
           ranges,
           fcontext,
           ppi_db,
           threads = 1) {
    
    # we set the OMP/BLAS number of threads to 1
    # this avoids issues we had in the glasso CV with multi-threading on cluster
    # also necessary for BDgraph
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    print("Fitting model using glasso.")
    glasso <- reg_net(data, priors, "glasso", threads = threads)
    
    print("Fitting bdgraph using priors.")
    bdgraph <- reg_net(data, priors, "bdgraph", threads = threads)
    
    print("Fitting model using iRafNet.")
    irafnet <- reg_net(data, priors, "irafnet", threads = threads)
    
    # --------------------------------------------------------------------------
    print("Add custom annotations for the graphs.")
    # --------------------------------------------------------------------------
    bdgraph$graph <- annotate.graph(bdgraph$graph,
                                    ranges, ppi_db, fcontext)
    
    irafnet$graph <- annotate.graph(irafnet$graph,
                                    ranges, ppi_db, fcontext)
    
    glasso$graph <- annotate.graph(glasso$graph,
                                   ranges, ppi_db, fcontext)
    
    # --------------------------------------------------------------------------
    print("Create result list.")
    # --------------------------------------------------------------------------
    result <- list(
      # bdgraph
      bdgraph_fit = bdgraph$fit,
      bdgraph = bdgraph$graph,
      # irafnet
      irafnet_fit = irafnet$fit,
      irafnet = irafnet$graph,
      # glasso
      glasso_fit = glasso$fit,
      glasso = glasso$graph
    )
    
    return(result)
  }



#' -----------------------------------------------------------------------------
#' Infers graphs for all available models and annotates them.
#' This method is used for testing subsets of models, i.e. we manually select
#' some models which should be run. At the moment these are genenet, glasso and
#' bdgraph without priors but using full start graph.
#'
#' @param data data matrix (n x p) to be used for fitting the models
#' @param priors the p x p prior matrix
#' @param ranges the ranges collection underlying the data
#' @param fcontext the tfbs context file
#' @param ppi_db the ppi db (graph) underlying the data
#' @param threads number of threads to be used. Default: 1
#'
#' -----------------------------------------------------------------------------
infer_all_graphs_subset <-
  function(data,
           priors,
           ranges,
           fcontext,
           ppi_db,
           threads = 1) {
    require(tidyverse)
    
    # we set the OMP/BLAS number of threads to 1
    # this avoids issues we had in the glasso CV with multi-threading on cluster
    # also necessary for BDgraph
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    print("Fitting model using GeneNet.")
    genenet <- reg_net(data, NULL, "genenet", threads = threads)
    
    # naive network, i.e. simply check corr > 0.5
    naive <- cor(data)
    naive[naive > 0.5] <- 1
    naive[naive != 1] <- 0
    naive_graph <- graphNEL(colnames(data))
    naive_graph <- with(
      melt(naive) %>% filter(value == 1),
      addEdge(as.character(Var1),
              as.character(Var2), naive_graph)
    )
    
    nlist <- list(naive = naive, graph = naive_graph)
    
    # --------------------------------------------------------------------------
    print("Add custom annotations for the graphs.")
    # --------------------------------------------------------------------------
    genenet$graph <- annotate.graph(genenet$graph,
                                    ranges, ppi_db, fcontext)
    
    nlist$graph <- annotate.graph(nlist$graph,
                                  ranges, ppi_db, fcontext)
    
    # --------------------------------------------------------------------------
    print("Create result list.")
    # --------------------------------------------------------------------------
    result <- list(
      # genenet
      genenet_fit = genenet$fit,
      genenet = genenet$graph,
      naive_fit = nlist$naive,
      naive = nlist$graph
    )
    
    return(result)
  }


#' -----------------------------------------------------------------------------
#' Test inference methods on a small, arbitrary ground truth graph and its data
#'
#' @param true_graph The ground truth graph
#' @param data data simulated according to the ground truth
#' @param iteration Iteration number (for annotating result frame)
#' @param threads number of threads to be used. Default: 1
#'
#' -----------------------------------------------------------------------------
test_inference <-
  function(true_graph, data, iteration, threads = 1) {
    require(tidyverse)
    
    # we set the OMP/BLAS number of threads to 1
    # this avoids issues we had in the glasso CV with multi-threading on cluster
    # also necessary for BDgraph
    RhpcBLASctl::omp_set_num_threads(1)
    RhpcBLASctl::blas_set_num_threads(1)
    
    print("Fitting model using GeneNet.")
    genenet <- reg_net(data, NULL, "genenet", threads = threads)
    
    print("Fitting model using glasso, no priors.")
    glasso_no_priors <-
      reg_net(data, NULL, "glasso", threads = threads)
    
    print("Fitting bdgraph without priors.")
    bdgraph_no_priors <- reg_net(data,
                                 NULL,
                                 "bdgraph",
                                 use_gstart = F,
                                 threads = threads)
    
    # naive network, i.e. simply check corr > 0.5
    naive <- cor(data)
    naive[naive > 0.5] <- 1
    naive[naive != 1] <- 0
    naive_graph <- graphNEL(colnames(data))
    naive_graph <- with(
      melt(naive) %>% filter(value == 1),
      addEdge(as.character(Var1),
              as.character(Var2), naive_graph)
    )
    
    nlist <- list(naive = naive, graph = naive_graph)
    
    # --------------------------------------------------------------------------
    print("Create result list.")
    result <- list(
      # naive
      naive_fit = nlist$naive,
      naive = nlist$graph,
      # bdgraph no priors
      bdgraph_no_priors_fit = bdgraph_no_priors$fit,
      bdgraph_no_priors = bdgraph_no_priors$graph,
      # genenet
      genenet_fit = genenet$fit,
      genenet = genenet$graph,
      # glasso no priors
      glasso_no_priors_fit = glasso_no_priors$fit,
      glasso_no_priors = glasso_no_priors$graph
    )
    
    gs <- result[!grepl("_fit", names(result))]
    
    # --------------------------------------------------------------------------
    # use the bdgraph internal method to get spec/sens, f1 and MCC. Use the
    # original simulation object containing the ground truth graph
    print("Validating...")
    
    perf <- lapply(names(gs), function(g) {
      # we need the adjacency matrix for comparison
      perf <-
        t(BDgraph::compare(as(true_graph, "matrix"), as(gs[[g]], "matrix")))
      comparisons <- c("True", g)
      perf <- as.data.frame(perf)
      rownames(perf) <- comparisons
      perf <- perf[!grepl("True", rownames(perf)),]
      
      # annotate density for comparison
      ig <- igraph::igraph.from.graphNEL(gs[[g]])
      dens <- edge_density(ig)
      perf$density_model <- dens
      perf$comparison = g
      perf
    }) %>% bind_rows() # lapply over different graph models --------------------
    
    # remember for easy plotting
    perf <- mutate(
      perf,
      snp = "random",
      iteration = iteration,
      density_true =
        edge_density(igraph.from.graphNEL(true_graph))
    )
    perf
  }


#' -----------------------------------------------------------------------------
#' Combine two graph objects.
#'
#' Keeps only nodes and edges present in both graphs
#'
#' @param g1 graphNEL object
#' @param g2 graphNEL object
#'
#' @import graph
#' @import igraph
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
combine_graphs <- function(g1, g2) {
  require(graph)
  require(igraph)
  
  # get incidence matrix for both graphs for easier comparison
  g1_mat <- as(g1, "matrix")
  g2_mat <- as(g2, "matrix")
  
  # use only common nodes, also ensure ordering of columns/rows
  use <- intersect(colnames(g1_mat), colnames(g2_mat))
  g1_mat <- g1_mat[use, use]
  g2_mat <- g2_mat[use, use]
  
  # get the graph from the combined matrix (i.e. only were edges are present in
  # both cohorts), also remove zero-degree nodes
  inc_combined <- (g1_mat == 1 & g2_mat == 1)
  g_combined <-
    igraph::as_graphnel(graph_from_adjacency_matrix(inc_combined,
                                                    mode = "undirected"))
  g_combined
}
