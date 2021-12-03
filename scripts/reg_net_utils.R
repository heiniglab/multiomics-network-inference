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

