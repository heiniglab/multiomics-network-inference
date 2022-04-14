# helper to quickly get a matrix containing edge information for a graph
get_edge_matrix <- function(g) {
  em <- t(edgeMatrix(g))
  n <- nodes(g)
  
  df <- cbind.data.frame(n1=n[em[,1]], n2=n[em[,2]], stringsAsFactors=F)
  # generate edges names (each name with sorted node names, to get easier
  # comparisons later on)
  rownames(df) <- sapply(1:nrow(df), function(x) { 
    paste0(sort(c(df[x,1], df[x,2])), collapse="~")
  })
  df
}

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

# generates merged graphs
# current strategy: annotate edges in the 'replicated' graph
# which are novel to the 'original' graph as well as 
# annotate 'replicated' edges
merge_graph <- function(novel_graph, original_graph) {
  
  nov_edges <- get_edge_matrix(novel_graph)
  ori_edges <- get_edge_matrix(original_graph)
  
  # define edges novel in ggm graph
  novel <- nov_edges[!rownames(nov_edges) %in% rownames(ori_edges), ]
  
  # define replicated edges
  replicated <-
    nov_edges[rownames(nov_edges) %in% rownames(ori_edges), ]
  
  # add edgeData
  g <- novel_graph
  edgeDataDefaults(g , "novel") <- F
  edgeDataDefaults(g, "replicated") <- F
  
  edgeData(g, replicated$n1, replicated$n2, "replicated") <- T
  edgeData(g, novel$n1, novel$n2, "novel") <- T
  
  g
}

# plotting function to specifically plot the merged graphs (recognizing 'replicated'
# edges)
# overall simpler than the other graph plots we had so far, since
# we do not annotate other edge types
plot_graph <- function(g, pdf_file, dot_file, ranges) {
  
  # get some default colors to be used here
  # we select custom colors for nice looks (and to fit to previous plots for the nodes)
  cols <- sort(get_defaultcolors(n=8))[c(1,5,3,7,2,6)]
  
  # set of node sin the graph
  n <- graph::nodes(g)
  
  #get genes for coloring, only keep the ones in our graph
  snp_genes <- ranges$snp_genes$SYMBOL
  tfs <- ranges$tfs$SYMBOL
  trans_genes <- ranges$cpg_genes$SYMBOL
  snp_genes <- snp_genes[snp_genes %in% n]
  tfs <- tfs[tfs %in% n]
  trans_genes <- trans_genes[trans_genes %in% n]
  
  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14,
                          style="filled", fontname="helvetica"),
                graph=list(overlap="false", outputorder="edgesfirst",
                           labelloc="top", labeljust="right"))
  
  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[grep("^cg", n)] = "box"
  shape[grep("^rs", n)] = "box"
  
  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("^cg", n)] = 0.4
  
  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("^cg", n)] = 0.4
  
  label = n
  names(label) = n
  #  label[grep("^cg", n)] = ""
  
  style <- rep("filled", numNodes(g))
  names(style) <- n
  
  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep("^rs", n)] = cols[1]
  col[grep("^cg", n)] = cols[2]
  
  if(length(snp_genes) > 0) {
    col[snp_genes] = cols[6]
  }
  if(length(tfs) > 0) {
    col[tfs] = cols[5]
  }
  if(length(trans_genes) > 0) {
    col[trans_genes] = cols[3]
  }
  
  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  
  bordercol = rep("black", numNodes(g))
  names(bordercol) = n
  bordercol[grep("^rs", n)] = "#ffe30f"
  
  # save node attributes
  nAttrs = list(shape=shape, label=label, style=style,  width=width,
                height=height, penwidth=penwidth, fillcolor=col,
                color=bordercol)
  
  # default color for edges: black
  ecol = rep("black", numEdges(g))
  names(ecol) = edgeNames(g)
  
  # check edgeData and add to colors
  for(edge in names(ecol)){
    n1 <- strsplit(edge,"~")[[1]][1]
    n2 <- strsplit(edge,"~")[[1]][2]
    
    if(unlist(graph::edgeData(g, n1, n2, "novel"))){
      ecol[edge] <- cols[3]
    }
    if(unlist(graph::edgeData(g,n1,n2, "replicated"))){
      ecol[edge] <- cols[4]
    }
  }
  
  dir = rep("none", numEdges(g))
  names(dir) = edgeNames(g)
  
  # save edge attributes
  eAttrs = list(color=ecol, dir=dir)
  
  # save graph as pdf
  pdf(pdf_file)
  plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  
  # Add legend, start with empty plot
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  
  legend("left", legend = c("SNP", "CpG", "gene", "novel edge","replicated edge"),
         pch=16, bg="#ff0000", pt.cex=3, cex=1.5, bty='n',
         col = c(cols[1], cols[2], "#ffffff", cols[3], cols[4]), 
         title="graph legend")
  
  dev.off()
  
  # and as dot for manual layouting
  toDot(g, dot_file, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
}