#' -----------------------------------------------------------------------------
#' Generate tables we show in the manuscript
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed Nov 20 13:11:37 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(igraph)
library(graph)
library(Rgraphviz)
library(BDgraph)

# ------------------------------------------------------------------------------
print("Table comparing our graphs to the ones reported in the meQTL study")
# ------------------------------------------------------------------------------
loci <- c("rs9859077", "rs7783715", "rs730775")

# get a graphNEL object from a dot file
parse_dot <- function(dot_file) {
  ag <- agread(dot_file)
  print(ag)
  enames <- Rgraphviz::edgeNames(ag)
  e1 <- unlist(lapply(strsplit(enames, "~"), "[[", 1))
  e2 <- unlist(lapply(strsplit(enames, "~"), "[[", 2))
  
  nodes <- unique(c(e1,e2))
  
  g <- graphNEL(nodes)
  g <- addEdge(e1,e2,g)
  g
}

# load our graphs
g1 <- sapply(loci, function(l) {
  parse_dot(paste0("results/current/biogrid_stringent/graph_plots_tfa/",
                   l,
                   "_meqtl/glasso_combined.dot"))
})

# load the ones from the meQTL paper
g2 <- sapply(loci, function(l) {
  parse_dot(paste0("../meQTLs/results/current/rw_ggm_integration/",
                   l,
                   "/graph-final.dot"))
})

# compare graphs
tab <- lapply(loci, function(l) {
  g2m <- as(g2[[l]], "matrix")
  g1m <- as(g1[[l]], "matrix")
  nodes_g2 <- colnames(g2m)
  nodes_g1 <- colnames(g1m)
  total_edges_g2 <- graph::numEdges(g2[[l]])
  total_edges_g1 <- graph::numEdges(g1[[l]])
  total_nodes_g1_in_g2 <- sum(nodes_g1 %in% nodes_g2)
  total_nodes_g2_in_g1 <- sum(nodes_g2 %in% nodes_g1)
  node_intersection <- intersect(nodes_g2, nodes_g1)
  
  g2m <- g2m[node_intersection, node_intersection]
  g1m <- g1m[node_intersection, node_intersection]
  
  # gets MCC and edge overlap (as true positives)
  comp <- BDgraph::compare(g2m, g1m)
  
  # selected SNP gene in referene
  snp_genes_in_ref <- unlist(adj(g2[[l]], l))
  snp_genes_in_est <- unlist(adj(g1[[l]], l))
  snp_genes_recovered <- intersect(snp_genes_in_ref, snp_genes_in_est)
  
  # return some stats
  rbind(c(locus=l,
    nodes_in_reference=length(nodes_g2),
    nodes_in_estimate=length(nodes_g1),
    edges_in_reference=total_edges_g2,
    edges_in_estimate=total_edges_g1,
    node_overlap=total_nodes_g2_in_g1,
    edge_overlap=comp["true positive", "estimate1"],
    snp_genes_reference=length(snp_genes_in_ref),
    snp_genes_estimate=length(snp_genes_in_est),
    snp_genes_recovered=length(snp_genes_recovered),
    mcc=comp["MCC", "estimate1"]))
})
tab <- do.call(rbind, tab) %>% 
  as_tibble()

# ----------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
