
ifile <- snakemake@input[[1]]
load(ifile)
ofile <- snakemake@output[[1]]

tab <- c()
temp <- lapply(names(simulated_data), function(n) {
  sd <- simulated_data[[n]]
  gh <- sd$graph.hidden
  gh_num_edges <- graph::numEdges(gh)
  gh_num_nodes <- graph::numNodes(gh)
  go <- sd$graph.observed
  go_num_edges <- graph::numEdges(go)
  go_num_nodes <- graph::numNodes(go)
  fnr <- sd$fnr
  fpr <- sd$fpr
  snp <- sd$snp
  
  tab <<- rbind(tab, c(snp,
                       gh_num_edges, gh_num_nodes,
                       go_num_edges, go_num_nodes,
                       fpr, fnr))
  colnames(tab) <- c("snp", "gh_num_edges", "gh_num_nodes",
                     "go_num_edges", "go_num_nodes", "fpr", "fnr")
})

# write to output file
write.table(file=ofile, tab, 
            col.names = T, row.names = F, sep="\t", quote=F)