#'
#' Defines functions overwriting library functions from
#' the BDgraph package in order to get what we think are more informative
#' plots for the graph results.
#'
#' @author  Johann Hawe
#' 


#' Here we define our own summary function for bdgraphs in order
#' to be able to avoid the graph plotting for large graphs (i.e. we
#' just add a flag for the original method  on whether or not to plot the graph...)
#'
#' @param object The ggm fit which to plot the summary for
#' @param vis Flag whether to visualize results
#' @param plot.graph Flag whether to plot the igraph or not (default:F) in case
#' vis = TRUE
#'
#' @author Johann Hawe
#'
summary.bdgraph <-function (object, plot.graph=F, vis = TRUE, ...) 
{
  p_links = object$p_links
  p = nrow(object$last_graph)
  dimlab = colnames(object$last_graph)
  selected_g = matrix(0, p, p, dimnames = list(dimlab, dimlab))
  if (!is.null(object$graph_weights)) {
    sample_graphs = object$sample_graphs
    graph_weights = object$graph_weights
    max_gWeights = max(graph_weights)
    sum_gWeights = sum(graph_weights)
    max_prob_G = max_gWeights/sum_gWeights
    if (is.null(dimlab)) 
      dimlab <- as.character(1:p)
    vec_G <- c(rep(0, p * (p - 1)/2))
    indG_max <- sample_graphs[which(graph_weights == max_gWeights)]
    vec_G[which(unlist(strsplit(as.character(indG_max), "")) == 
                  1)] = 1
    selected_g[upper.tri(selected_g)] <- vec_G
  } else {
    selected_g[p_links > 0.5] = 1
    selected_g[p_links <= 0.5] = 0
  }
  if (vis) {
    G <- graph.adjacency(selected_g, mode = "undirected", 
                         diag = FALSE)
    if (!is.null(object$graph_weights)) {
      op = par(mfrow = c(2, 2), pty = "s", omi = c(0.3, 
                                                   0.3, 0.3, 0.3), mai = c(0.3, 0.3, 0.3, 0.3))
      subGraph = paste(c("Posterior probability = ", max_prob_G), 
                       collapse = "")
    }
    else {
      subGraph = "Selected graph with edge posterior probability = 0.5"
    }
    if (p < 20) 
      size = 15
    else size = 2
    if(plot.graph) {
      plot.igraph(G, layout = layout.circle, main = "Selected graph", 
                  sub = subGraph, vertex.color = "white", vertex.size = size, 
                  vertex.label.color = "black")
    }
    if (!is.null(object$graph_weights)) {
      plot(x = 1:length(graph_weights), y = graph_weights/sum_gWeights, 
           type = "h", main = "Posterior probability of graphs", 
           ylab = "Pr(graph|data)", xlab = "graph")
      abline(h = max_prob_G, col = "red")
      text(which(max_gWeights == graph_weights)[1], max_prob_G, 
           "Pr(selected graph|data)", col = "gray60", adj = c(0, 
                                                              +1))
      sizesample_graphs = sapply(sample_graphs, function(x) length(which(unlist(strsplit(as.character(x), 
                                                                                         "")) == 1)))
      xx <- unique(sizesample_graphs)
      weightsg <- vector()
      for (i in 1:length(xx)) weightsg[i] <- sum(graph_weights[which(sizesample_graphs == 
                                                                       xx[i])])
      plot(x = xx, y = weightsg/sum_gWeights, type = "h", 
           main = "Posterior probability of graphs size", 
           ylab = "Pr(graph size|data)", xlab = "Graph size")
      all_graphs = object$all_graphs
      sizeall_graphs = sizesample_graphs[all_graphs]
      plot(x = 1:length(all_graphs), sizeall_graphs, type = "l", 
           main = "Trace of graph size", ylab = "Graph size", 
           xlab = "Iteration")
      abline(h = sum(selected_g), col = "red")
      par(op)
    }
  }
  if (!is.null(object$graph_weights)) {
    pvec <- 0 * vec_G
    for (i in 1:length(sample_graphs)) {
      which_edge <- which(unlist(strsplit(as.character(sample_graphs[i]), 
                                          "")) == 1)
      pvec[which_edge] <- pvec[which_edge] + graph_weights[i]
    }
    p_links <- 0 * selected_g
    p_links[upper.tri(p_links)] <- pvec/sum_gWeights
  }
  K_hat = object$K_hat
  if (is.null(K_hat)) 
    return(list(selected_g = Matrix(selected_g, sparse = TRUE), 
                p_links = Matrix(p_links, sparse = TRUE)))
  else return(list(selected_g = Matrix(selected_g, sparse = TRUE), 
                   p_links = Matrix(p_links, sparse = TRUE), K_hat = K_hat))
}
