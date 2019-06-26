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

# define the available models for network inference
reg_net.models <- function() {
  return(c("genenet", "bdgraph", "irafnet", "glasso", "genie3" ,"custom"))
}

# ------------------------------------------------------------------------------
#' Main method for inferring regulatory networks.
#'
#' @param data The data matrix (n x p) from which to infer the network
#' @param priors The matrix of priors (p x p) which to use. Needs to be set to
#' NULL explicitely. In case of bdgraph model, these are used to define a start
#' graph for the algorithm.
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
#' @param irafnet.fdr optional. FDR cutoff for edges in final iRafNet graph.
#' Default: 0.05
#'
#' @return Returns a list containing bot the final model fit as well as the
#' graphNEL object extracted from that model fit.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
# ------------------------------------------------------------------------------
reg_net <- function(data, priors, model, threads=1,
                    use_gstart=T, gstart=NULL, iter=10000, burnin=5000,
                    ntrees=1000, mtry=round(sqrt(ncol(data)-1)), npermut=5,
                    irafnet.fdr=0.05, glasso.lambda=1) {

  # load inference methods
  suppressPackageStartupMessages(library(GeneNet))
  suppressPackageStartupMessages(library(BDgraph))
  suppressPackageStartupMessages(library(iRafNet))
  suppressPackageStartupMessages(library(glasso))
  suppressPackageStartupMessages(library(GENIE3))

  # get available models
  ms <- reg_net.models()
  if(!(model %in% ms)) stop(paste0("Model not supported: ", model))

  # for genenet and iRafNet, remove NAs
  data_no_nas <- data[,apply(data, 2, function(x) !anyNA(x))]

  if(!is.null(priors)) {
    # we possibly lost some data when filtering for NAs, so we adjust priors too
    priors_no_nas <- priors[rownames(priors) %in% colnames(data_no_nas),
                            colnames(priors) %in% colnames(data_no_nas)]
  } else {
    priors_no_nas <- NULL
  }
  # check which model to build
  if("genenet" %in% model){
    pcors <- ggm.estimate.pcor(data_no_nas)
    fit <- network.test.edges(pcors, plot = F)
    fit$node1 <- colnames(data_no_nas)[fit$node1]
    fit$node2 <- colnames(data_no_nas)[fit$node2]
    class(fit) <- c(class(fit), "genenet")
  } else if("bdgraph" %in% model) {
    # do we have priors?
    if(is.null(priors)) {
      # build without priors, nevertheless, we could have a custom start graph
      if(use_gstart) {
        if(is.null(gstart)) {
          stop("Start graph must be provided if no priors are given.")
        }
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin, g.start=gstart,
                       save = T, cores=threads)
      } else {
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       save = T, cores=threads)
      }
    } else {
      # we have priors, check whether to use the prior based start graph
      if(use_gstart) {
        gstart <- get_gstart_from_priors(priors)
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       g.prior = priors, g.start = gstart,
                       save = T, cores=threads)
      } else {
        # use priors, but no specific start graph
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       g.prior = priors, g.start="empty",
                       save = T, cores=threads)
      }
    }

    # plot convergence info for any case
    ggm_summary <- summary(fit)

  } else if("irafnet" %in% model) {
    irn_out <- iRafNet(data_no_nas, priors_no_nas, ntrees, mtry, colnames(data_no_nas),
                       threads=threads)
    irn_perm_out <- Run_permutation(data_no_nas, priors_no_nas,
                                    ntrees, mtry, colnames(data_no_nas),
                                    npermut, threads=threads)
    fit <- iRafNet_network(irn_out, irn_perm_out, TH = irafnet.fdr)
    fit <- list(fit=fit, nodes=colnames(data_no_nas))
    class(fit) <- c(class(fit), "irafnet")
  } else if("custom" %in% model) {
    stop("Sorry, custom model is not yet implemented.")
  } else if("glasso" %in% model) {
    if(!is.null(priors)) {
      # for now we simply call using 1-priors as penalization
      gl_out <- glasso_cv(data, priors, threads = threads)
    } else {
      gl_out <- glasso_cv(data, threads = threads)
    }
    class(gl_out) <- c(class(gl_out), "glasso")
    fit <- gl_out
  } else if("genie3" %in% model) {
    fit <- genie3(data_no_nas, threads = threads)
    class(fit) <- c(class(fit), "genie3")
  }

  # now get the graph object
  g <- graph_from_fit(fit, annotate = F)

  return(list(graph=g, fit=fit))
}

# ------------------------------------------------------------------------------
#' Creates a graphNEL object from a given bdgraph result for a defined cutoff
#'
#' @param ggm.fit The bdgraph ggm fit
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
                           ranges = NULL,
                           ppi_db=NULL,
                           annotate=T){

  if(annotate & is.null(ranges)) {
    stop("Ranges must not be null for annotating graphs!")
  }

  suppressPackageStartupMessages(library(BDgraph))
  suppressPackageStartupMessages(library(graph))
  suppressPackageStartupMessages(library(igraph))
  require(reshape2)
  #  suppressPackageStartupMessages(require(tidyverse))

  # get the graph instance from the ggm fits
  if (inherits(ggm.fit, "bdgraph")) {
    g.adj <- BDgraph::select(ggm.fit, cut = 0.9)
    g <-
      as_graphnel(graph.adjacency(g.adj, mode = "undirected", diag = F))
  } else if (inherits(ggm.fit, "irafnet")) {
    g <- graphNEL(ggm.fit$nodes, edgemode = "undirected")
    fit <- ggm.fit$fit
    g <- addEdge(fit$gene1, fit$gene2, g)
  } else if (inherits(ggm.fit, "genenet")) {
    n <- unique(c(ggm.fit$node1, ggm.fit$node2))
    net <- extract.network(ggm.fit,
                           cutoff.ggm = 0.8)
    g <- graphNEL(n,
                  edgemode = "undirected")
    g <- addEdge(net$node1, net$node2, g)
  } else if(inherits(ggm.fit, "glasso")) {
    pm <- ggm.fit$wi
    g <- graphNEL(colnames(pm),
                  edgemode="undirected")
    # create edge matrix
    temp <- melt(pm)

    # convert back to characters for graphNEL
    temp$Var1 <- as.character(temp$Var1)
    temp$Var2 <- as.character(temp$Var2)

    # define edges and remove self-edges
    temp <- subset(temp, Var1 != Var2 & value != 0)

    if(nrow(temp) > 0) {
      g <- addEdge(temp$Var1, temp$Var2, g)
    }
  } else if(inherits(ggm.fit, "genie3")) {
    g <- get_genie3_graph(ggm.fit$nodes,
                          ggm.fit$model,
                          ggm.fit$best_weight)
  }

  if(annotate) {
    # set node and edge attributes
    g <- annotate.graph(g, ranges, ppi_db)
  }
  return(g)
}

# ------------------------------------------------------------------------------
#' Fits the glasso models for a selection of lambdas and identifies the one
#' no priors fit which is closest (edge number wise) to the prior model with
#' no scaling
#'
#' DEPRECATED: use glasso_cv instead!
# ------------------------------------------------------------------------------
glasso_screen <- function(data, priors, threads, ranges, ppi_db, fcontext) {
  warning("Method glasso_screen is deprecated. Use glasso_cv instead!")

  require(doParallel)

  lambdas <- seq(0.1,1,by=0.1)
  res <- foreach(l = lambdas, .export=c("reg_net", "reg_net.models", "annotate.graph", "graph_from_fit",
                                        "filter.edge.matrix", "get_tfbs_context"),
                 .packages=c("reshape2", "glasso", "graph")) %dopar% {
                   glasso <- reg_net(data, priors, "glasso", threads=threads, glasso.lambda=l)
                   glasso$graph <- annotate.graph(glasso$graph, ranges, ppi_db, fcontext)

                   # here we use the average penalty arising from the scaled prior matrix
                   ut <- l * (1-priors[upper.tri(priors)])
                   lambda <- sum(ut)/length(ut)
                   glasso_no_priors <- reg_net(data, NULL, "glasso", glasso.lambda=lambda, threads=threads)
                   glasso_no_priors$graph <- annotate.graph(glasso_no_priors$graph, ranges, ppi_db, fcontext)


                   glassos <- list()
                   n <- paste0("glasso_lambda", l)
                   n2 <- paste0(n, "_no_priors")
                   glassos[[n]] <- glasso
                   glassos[[n2]] <- glasso_no_priors
                   return(glassos)
                 }
  names(res) <- paste0("lambda", lambdas)

  # we directly get the two main glasso results, i.e. the one with no modulation
  # to the priors and the one with no priors which has similar number of edges
  # NOTE: for additional evaluation possibilities, we also save all of the trained
  # lasso models

  gl <- res[["lambda1"]]$glasso_lambda1
  ne <- numEdges(gl$graph)
  closest_graph <- NA
  edge_diff <- NA

  for(l in names(res)) {
    # select corresponding non_prior_graph
    re <- res[[l]]
    gnp <- re[[which(grepl(".*_no_priors", names(re)))]]
    ed <- abs(ne - numEdges(gnp$graph))
    if(is.na(edge_diff)) {
      edge_diff <- ed
      closest_graph <- gnp
    } else {
      if(ed < edge_diff) {
        edge_diff <- ed
        closest_graph <- gnp
      }
    }
  }
  glnp <- closest_graph

  glasso_all <- res
  return(list(glasso_all=glasso_all, glasso_no_priors = glnp, glasso = gl))
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
get_graph_score <- function(g, sentinel, ranges, density=NULL) {
  require(igraph)

  seed <- ranges$seed

  # get graph score
  score <- 0

  # keep only cluster with sentinel
  cl <- clusters(g)
  if(!sentinel %in% names(cl$membership)) return(0)

  cn <- cl$membership[sentinel]
  subnodes <- names(cl$membership[cl$membership == cn])
  if(length(subnodes) < 2) return(0)

  g <- induced_subgraph(g, subnodes)

  # check whether we have any SNP gene.
  # if not, the score will be 0
  # otherwise proceed with check
  adj_sent <- igraph::neighbors(g, sentinel)$name

  # get number of trans genes in our cluster
  v <- V(g)$name
  if(seed == "meqtl") {
    tg_in_v <- intersect(names(ranges$cpgs), v)
  } else {
    tg_in_v <- intersect(ranges$trans_genes$SYMBOL, v)
  }

  if(seed %in% c("meqtl", "eqtlgen", "eqtl")) {

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
    Y_trans <- 1/length(tg_in_v)

    for(trans in tg_in_v) {
      # remove all other trans-genes to get only independently reachable TGs
      sg <- delete.vertices(g, setdiff(tg_in_v, trans))

      paths <-
        suppressWarnings(get.shortest.paths(sg, sentinel, trans))$vpath[[1]]

      # we need at least one TF or Spath gene on the path
      pgenes <- unique(c(ranges$tfs$SYMBOL, ranges$spath$SYMBOL))

      if(length(paths) > 1 && any(pgenes %in% paths$name)) {
        score <- score + Y_trans
      }
    }
  }
  # if graph density is available, we use it to adjust the score
  if(!is.null(density)) {
    score <- score * (1-density)
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
  c2 <- p*log(2*pi)
  logdet <- determinant(theta, logarithm=T)$modulus
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
  return( log(n) * k - 2 * ll)
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
glasso_cv <- function(data, priors = NULL, k=5,
                      rholist=seq(0.01,1,by=0.005),
                      threads = 1) {
  require(glasso)
  require(cvTools)
  require(parallel)

  n <- nrow(data)
  folds <- cvFolds(n, k, type="random")

  # get the ll progression for all CV folds
  loglikes <- lapply(1:k, function(ki) {
    S_train <- cov(data[folds$which!=ki,], use="na.or.complete")
    S_test  <- cov(data[folds$which==ki,], use="na.or.complete")

    # get the number of variables
    p <- ncol(S_test)
    n <- nrow(S_test)

    # check all rhos
    # NOTE: we can't use glassopath since we want to include the prior information
    unlist(mclapply(rholist, function(rho) {
      if(!is.null(priors)) {
        gl <- glasso(S_train, rho = rho*(1-priors),
                     penalize.diagonal = F)
      } else {
        gl <- glasso(S_train, rho = rho,
                     penalize.diagonal = F)
      }
      rownames(gl$wi) <- colnames(gl$wi) <- colnames(data)
      class(gl) <- c(class(gl), "glasso")

      # get log likelihood, number of edges and BIC
      ll <- loglik(gl$wi, S_test, n, p)

      g <- graph_from_fit(gl, annotate=F)
      ne <- numEdges(g)

      bic <- bic(ll=ll, n=n, k = ne)

      c(ll, bic, ne)

    }, mc.cores=threads))
  })
  loglikes <- do.call(rbind, loglikes)
  colnames(loglikes) <- paste0(c("rho_rho=", "bic_rho=", "edgeNum_rho="),
                               rep(rholist, each=3))
  rownames(loglikes) <- paste0("k=", 1:k)

  # get the 'best' rho according to the BIC
  bic_sub <- colMeans(loglikes[, grepl("bic", colnames(loglikes))])
  ind <- which.min(bic_sub)

  # retrieve rho from the name
  rho <- as.numeric(strsplit(names(ind), "=")[[1]][2])

  # retrain the full model
  S <- cov(data, use="na.or.complete")
  if(!is.null(priors)) {
    model <- glasso(S, rho*(1-priors), penalize.diagonal = F)
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
#' @param nodes All nodes used in the calculation of the model
#' @param model The trained GENIE3 model
#' @param threshold The treshold to be applied in GENIE3's \code{getLinkList}
#'
#' -----------------------------------------------------------------------------
get_genie3_graph <- function(nodes, model, threshold) {
  require(graph)
  require(GENIE3)

  edges <- getLinkList(model, threshold = threshold)
  if(nrow(edges) > 0) {
    g <- graph::graphNEL(nodes,
                  edgemode="undirected")
    g <- graph::addEdge(as.character(edges$regulatoryGene),
                        as.character(edges$targetGene), g)
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
#' for a specific threshold to a powerlow and then selecting the threshold for
#' which the likelihood of the fitted model is maximal
#'
#' @param data The data matrix (n x p)
#' @param threads The number of threads to be used
#'
#' -----------------------------------------------------------------------------
genie3 <- function(data, threads=1) {
  require(GENIE3)
  require(igraph)
  require(parallel)

  print(paste0("Using ", threads, " threads."))

  # we expect a n x p matrix, genie3 needs a p x n matrix
  model <- GENIE3(t(data), nCores = threads)

  all_link_weights <- sort(unique(getLinkList(model, threshold = 0)$weight))
  print("Summary of link weights:")
  print(summary(all_link_weights))

  lw <- length(all_link_weights)
  # get a subset to check if needed
  if(lw > 500) {
    all_link_weights <- quantile(all_link_weights, seq(0,1,by = 0.005))
  }

  print(paste0("Checking ", length(all_link_weights), " weights."))
  print("Summary of link weights:")
  print(summary(all_link_weights))

  # we can use igraph packages 'fit_power_law' to fit the degree distribution
  # to a power law distr for each weight cutoff
  fits <- mclapply(all_link_weights, function(weight) {
    g <- get_genie3_graph(colnames(data), model, weight)
    ds <- graph::degree(g)
    # fit pwer law
    fit <- igraph::fit_power_law(ds+1)
    # return ll of fitted distr
    c(ll=fit$logLik, ks_p=fit$KS.p)
  }, mc.cores=threads)

  lls <- unlist(lapply(fits, "[[", "ll"))
  names(lls) <- paste0("weight=", all_link_weights)
  ks_p <- unlist(lapply(fits, "[[", "ks_p"))
  names(ks_p) <- paste0("weight=", all_link_weights)

  # get best weight (highest KS_p and heighest weight for that)
  best_weights <- all_link_weights[lls == max(lls)]
  best_weight <- best_weights[which.max(best_weights)]

  return(list(model=model, nodes=colnames(data),
              best_weight=best_weight, logLiks=lls, KS_ps=ks_p))
}
