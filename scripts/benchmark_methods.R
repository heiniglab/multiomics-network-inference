#' -----------------------------------------------------------------------------
#' Methods script for the benchmarking.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Apr  8 07:46:06 2021
#' -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Creates a simulated data and prior matrix to be used for the benchmarking
#'
#' @param number_of_nodes The number of nodes to simulate  
#' @param number_of_samples The number of samples to simulate
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
simulate_data <- function(number_of_nodes, number_of_samples) {
  require(BDgraph)
  
  node_names <- paste0("N", 1:number_of_nodes)
  
  data <- BDgraph::bdgraph.sim(p=number_of_nodes, 
                               n = number_of_samples, 
                               graph = "scale-free")$data
  colnames(data) <- node_names
  
  
  priors <- create_random_prior_matrix(data)
  colnames(priors) <- rownames(priors) <- node_names
  
  return(list(data = data, priors = priors))
}

# ------------------------------------------------------------------------------
#' Method DESC
#'
#' @param 
#' @param
#' @param
#' @param
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
create_random_prior_matrix <- function(data) {
  
  number_of_nodes <- ncol(data)
  
  PSEUDO_PRIOR <- 1e-7
  priors <- matrix(PSEUDO_PRIOR, nrow=number_of_nodes, ncol=number_of_nodes)
  
  # set random prior values
  number_of_priors <- floor(runif(1) * number_of_nodes)
  prior_indices <- sample(which(upper.tri(priors)), number_of_priors)
  
  priors[prior_indices] <- runif(number_of_priors, PSEUDO_PRIOR, 1 - PSEUDO_PRIOR)
  
  # make symmetric
  priors[lower.tri(priors)] <- t(priors)[lower.tri(priors)]
  
  return(priors)
}
