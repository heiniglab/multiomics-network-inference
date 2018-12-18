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
  return(c("genenet", "bdgraph", "irafnet", "custom"))
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
                    use_gstart=T, gstart=NULL, iter=10000, burnin=2500,
                    ntrees=1000, mtry=round(sqrt(ncol(data)-1)), npermut=5,
                    irafnet.fdr=0.05) {

  # get available models
  ms <- reg_net.models()
  if(!(model %in% ms)) stop(paste0("Model not supported: ", model))

  # check which model to build
  if("genenet" %in% model){
    # for genenet, remove NAs
    gn_data <- data[,apply(data, 2, function(x) !anyNA(x))]
    pcors <- ggm.estimate.pcor(gn_data)
    fit <- network.test.edges(pcors, plot = F)
    fit$node1 <- colnames(gn_data)[fit$node1]
    fit$node2 <- colnames(gn_data)[fit$node2]
    class(fit) <- c(class(fit), "genenet")
  } else if("bdgraph" %in% model) {
    # do we have priors?
    if(is.null(priors)) {
      # build without priors, nevertheless, we could have a custom start graph
      if(use_gstart) {
        if(is.null(gstart)) {
          stop("Start graph must be provided if now priors are given.")
        }
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin, g.start=gstart,
                       save.all=T, cores=threads)
      } else {
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       save.all=T, cores=threads)
      }
    } else {
      # we have priors, check whether to use the prior based start graph
      if(use_gstart) {
        gstart <- get_gstart_from_priors(priors)
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       g.prior = priors, g.start = gstart,
                       save.all=T, cores=threads)
      } else {
        # use priors, but no specific start graph
        fit <- bdgraph(data, method = "gcgm",
                       iter = iter, burnin = burnin,
                       g.prior = priors, g.start="empty",
                       save.all=T, cores=threads)
      }
    }
  } else if("irafnet" %in% model) {
    irn_out <- iRafNet(data, priors, ntrees, mtry, colnames(data),
                       threads=threads)
    irn_perm_out <- Run_permutation(data, priors,
                                    ntrees, mtry, colnames(data),
                                    npermut, threads=threads)
    fit <- iRafNet_network(irn_out, irn_perm_out, TH = irafnet.fdr)
    class(fit) <- c(class(fit), "irafnet")
  } else if("custom" %in% model) {
    stop("Sorry, custom model is not yet implemented.")
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
#' @param fcontext
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
                           string_db=NULL,
                           fcontext=NULL,
                           annotate=T){

  if(annotate & (is.null(fcontext) | is.null(ranges))) {
    stop("Chip-seq context and ranges must not be null for annotating graphs!")
  }

  suppressPackageStartupMessages(library(BDgraph))
  suppressPackageStartupMessages(library(graph))
  suppressPackageStartupMessages(library(igraph))

  # get the graph instance from the ggm fit
  cutoff <- 0.95
  if (inherits(ggm.fit, "bdgraph")) {
    g.adj <- BDgraph::select(ggm.fit, cut = cutoff)
    g <-
      as_graphnel(graph.adjacency(g.adj, mode = "undirected", diag = F))
  } else if (inherits(ggm.fit, "irafnet")) {
    g <- graphNEL(unique(unlist(ggm.fit)), edgemode = "undirected")
    g <- addEdge(ggm.fit$gene1, ggm.fit$gene2, g)
  } else if (inherits(ggm.fit, "genenet")) {
    net <- extract.network(ggm.fit,
                           cutoff.ggm = cutoff)
    g <- graphNEL(unique(c(net$node1, net$node2)),
                  edgemode = "undirected")
    g <- addEdge(net$node1, net$node2, g)

  }

  if(annotate) {
    # set node and edge attributes
    g <- annotate.graph(g, ranges, string_db, fcontext)
  }
  return(g)
}
