#' Get human gene annotation with symbols.
#'
#' This will get the gene annotation annotated with respective gene SYMBOL from UCSC for the hg19 build
#' TODO: allow other builds as well, currently this method is not very flexible...
#'
#' @param drop.nas Flags whether to drop any annotations where the SYMBOL turns out to be NA
#' @param verbose Verbosity flag (default: TRUE)
#'
#' @return GRanges object, where names(obj) are entrez-ids, and obj$SYMBOLS contains respective gene symbols
#'
get.gene.annotation <- function(drop.nas=TRUE, verbose=TRUE) {
  library(GenomicFeatures)
  library(GenomicRanges)
  library(org.Hs.eg.db)
  
  # TODO we'd like to use the faster way as indicated below. However this yields
  # roughly 400 less genes as compared to downloading the current version from
  # UCSC directly -> we use the list with more genes for now
  # (yes, it is also the difference between genes/individual transcripts, but some
  # gene symbols gotten via the current approach are not found when using the easier
  # one...)
  #  require(TxDb.Hsapiens.UCSC.hg19.knownGene
  #  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  txdb <- makeTxDbFromUCSC(genome="hg19",tablename="knownGene")
  
  tscripts = transcriptsBy(txdb, by="gene")

  # names(gene.ranges) gives entrezid
  gene.ranges <- unlist(range(tscripts))

  gene.ranges$SYMBOL <- unlist(mget(names(gene.ranges), org.Hs.egSYMBOL, ifnotfound=NA))
  if(drop.nas){
    to.keep <- which(!is.na(gene.ranges$SYMBOL));
    if(verbose) cat("Removing", length(gene.ranges)-length(to.keep), "genes since SYMBOL = NA.\n");
    gene.ranges <- gene.ranges[to.keep]
  }
  return(gene.ranges)
}

#'
#' Define a convenience function to get linear model p-values
#'
#' @param modelobject The linear model for which to get the p-value
#'
#' @return The p-value to the given linear model
#'
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

#' Gets a single snp range from the gtex snp database
#' 
#' @author Johann Hawe
#' 
get.snp.range <- function(snp){
  library(data.table)
  snps <- fread(paste0("results/current/gtex-snp-locations.txt"))
  snps <- snps[RS_ID_dbSNP142_CHG37p13 %in% snp,]
  if(nrow(snps)<1){
    warning("SNP not found on gtex db.\n") 
    return(NULL)
  } else {
    r <- with(snps, 
              GRanges(paste0("chr", Chr), 
                      IRanges(as.integer(Pos),width=1)))
    names(r) <- snps$RS_ID_dbSNP142_CHG37p13
    return(r)
  }
}

#' Gets the chip-seq context for a specific set of cpgs
#' 
#' @param cpgs ID of cpgs for which to get the context
#'
#' @author Johann Hawe
#' 
get.chipseq.context <- function(cpgs){
  load("data/cpgs_with_chipseq_context_100.RData")
  tfbs.ann <- tfbs.ann[rownames(tfbs.ann) %in% cpgs,,drop=F]
  
  return(tfbs.ann[cpgs,,drop=F])
}

#' Creates a graphNEL object from a given bdgraph result for a defined cutoff
#' 
#' @param ggm.fit The bdgraph ggm fit
#' @param ranges The ranges of the entities used for the graph fit
#' 
#' @value graph-nel object created from the bdgraph result
#'
#' @author Johann Hawe
#'
graph.from.fit <- function(ggm.fit, ranges){
  library(BDgraph)
  library(graph)
  library(igraph)
  
  # get the graph instance from the ggm fit
  cutoff <- 0.8
  g.adj <- BDgraph::select(ggm.fit, cut = cutoff)
  g <- as_graphnel(graph.adjacency(g.adj, mode="undirected", diag=F));
  gn <- nodes(g)
  
  nodeDataDefaults(g,"cpg") <- F
  nodeData(g, gn, "cpg") <- grepl("^cg", gn)
  
  nodeDataDefaults(g,"snp") <- F
  nodeData(g, gn, "snp") <- grepl("^rs", gn)
  
  nodeDataDefaults(g,"snp.gene") <- F
  nodeData(g, gn, "snp.gene") <- gn %in% ranges$snp.genes$SYMBOL
  
  nodeDataDefaults(g,"cpg.gene") <- F
  nodeData(g, gn, "cpg.gene") <- gn %in% ranges$cpg.genes$SYMBOL
  
  nodeDataDefaults(g, "tf") <- F
  nodeData(g, gn, "tf") <- gn %in% ranges$enriched.tfs$SYMBOL
  
  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% ranges$shortestpath.genes$SYMBOL
  
  # edge data information
  edgeDataDefaults(g, "isChipSeq") <- FALSE
  edgeDataDefaults(g, "isPPI") <- FALSE 
  
  # we will also set tf2cpg priors at this point, so get all TFs
  tfs <- ranges$enriched.tfs$SYMBOL
  # also get the chipseq context for our cpgs
  context <- get.chipseq.context(names(ranges$cpgs))
  
  em <- matrix(ncol=2,nrow=0)
  # for all cpgs
  for(c in rownames(context)){
    for(tf in tfs) {
      # might be that the TF was measured in more than one cell line
      if(any(context[c,grepl(tf, colnames(context))])) {
       em <- rbind(em,c(c,tf))
      }
    }
  }
  if(nrow(em) > 0){
    em <- filter.edge.matrix(g,em)

    edgeData(g,em[,1], em[,2],"isChipSeq") <- T
  }
  # ppi edgedata
  load.string.db()
   # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), gn), STRING.DB)
  edges <- t(edgeMatrix(STRING.SUB))
  edges <- cbind(gn[edges[,1]], gn[edges[,2]])
  edges <- filter.edge.matrix(g,edges)
  if(nrow(edges) > 0){
    edgeData(g,edges[,1], edges[,2],"isPPI") <- T
  }
  return(g)
}

#' Plot a GGM result graph
#'
#' Plots the a built graph (estimated from the sentinel data) using the 
#' twopi-visualization. Can be used to retrieve only the plot graph/node/edge
#' attributes when using pdf.out=NULL, dot.out=NULL and plot.on.device=F
#'
#' @param graph: Graph to be plotted (graphNEL)
#'
#' @param id The Id of the sentinel snp
#' @param ranges The list of sentinel associated ranges (we need snp.genes 
#' and cpg.genes, enriched.tfs etc.)
#' @param dot.out File to which to write the graph in dot format
#' 
#' @return List of plot graph attributes and the underlying graph structure
#' 
#' @author Johann Hawe
#'
plot.ggm <- function(g, id, dot.out=NULL){
  library(graph)
  library(Rgraphviz)
  library(GenomicRanges)
  
  # remove any unconnected nodes (primarily for bdgraph result, since
  # such nodes are already removed for genenet)
  if(any(graph::degree(g) == 0)){
    g <- removeNode(names(which(graph::degree(g) == 0)), g)
  }
  
  # add sentinel to network if its is not in there yet or it has been removed
  if(!(id %in% nodes(g))) {
    g <- graph::addNode(c(id), g)
  }
  
  n <- nodes(g)
  
  # get trans and cpg gene symbols
  snp.genes <- n[unlist(nodeData(g,n,"snp.gene"))]
  cpg.genes <- n[unlist(nodeData(g,n,"cpg.gene"))]
  tfs <- n[unlist(nodeData(g,n,"tf"))]
  
  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14, 
                          style="filled", fontname="helvetica"), 
                graph=list(overlap="false", root=id, outputorder="edgesfirst"))

  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[grep("^cg", n)] = "box"
  shape[grep("^rs", n)] = "box"

  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("cg", n)] = 0.4

  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("cg", n)] = 0.4

  label = n
  names(label) = n
  label[grep("cg", n)] = ""

  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep("^rs", n)] = "#fab4ad";
  col[grep("^cg", n)] = "#e4d7bc";
  if(!is.null(tfs)){
    col[tfs] = "green"
  }
  
  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  penwidth[snp.genes] = 3
  penwidth[cpg.genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3  
  }
 
  bordercol = rep("black", numNodes(g));
  names(bordercol) = n;
  bordercol[cpg.genes] = "#e4d7bc";
  bordercol[id] = "#fab4ad";

  nAttrs = list(shape=shape, label=label, width=width, 
                height=height, penwidth=penwidth, fillcolor=col, 
                color=bordercol)

  # default color for edges: black
  ecol = rep("black", numEdges(g))
  names(ecol) = edgeNames(g)
  for(i in snp.genes) {
    # color any edge from a SNP to one of our snp genes red
    ecol[grepl("^rs|~rs", names(ecol)) & grepl(i, names(ecol))] = "#b3cde2"
  }

  # set also color for cpgs
  for(cg in cpg.genes){
    # color any edge from a cpg to one of its cpg genes blue (proximity edges)
    ecol[grepl("^cg|~cg", names(ecol)) & grepl(cg, names(ecol))] = "#b3cde2"
  }

  # check edgeData and add to colors
  for(edge in names(ecol)){
    n1 <- strsplit(edge,"~")[[1]][1]
    n2 <- strsplit(edge,"~")[[1]][2]

    if(unlist(graph::edgeData(g,n1,n2, "isPPI"))){
      ecol[edge] <- "#decae3"
    }
    if(unlist(graph::edgeData(g,n1,n2, "isChipSeq"))){
      ecol[edge] <- "#ccebc5"
    }
  }
  
  dir = rep("none", numEdges(g))
  names(dir) = edgeNames(g)

  eAttrs = list(color=ecol, dir=dir)
 
  if(numEdges(g)>500){
    warning("Skipping plotting on device due to large amount of edges")
  } else{
    plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  if(!is.null(dot.out)){
    # output the dot-information
    #    toDot(graph, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs, subGList=sgs)
    toDot(g, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  # return the list of created plotattributes and the possibly modified graph 
  # object
  return(list(graph=g, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs))
}

#' Method to quickly filter an edge matrix for only those edges, which are within
#' a specified graph
#' 
#' @author Johann Hawe
#' 
#' @date 2017/03/13
#' 
filter.edge.matrix <- function(g, em){
  e <- graph::edges(g)
  out <- matrix(ncol=2,nrow=0)
  for(i in 1:nrow(em)){
    e1 <- em[i,1]
    e2 <- em[i,2]
    if(e1 %in% names(e)){
      if(e2 %in% e[[e1]]){
        out <- rbind(out,c(e1,e2))
      }
    }
  }
  return(out)
}

#' Get  STRING interactions (only experimental and db supported)
#'
#' Loads the STRING db saved on disc to a graph R object.
#' Uses a caching mechanism for faster loading of the graph information
#'
#' @return nothing
#'
load.string.db <- function(reload=F) {
  library(data.table)
  library(igraph)
  library(graph)
  
  cat("Loading string db\n")
 
  fcache <- "results/string.v9.expr.RData"
  if(reload || !file.exists(fcache)){
  # load db anew
  string.all <- fread(paste0("data/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt"),
                        data.table=F, header=T, stringsAsFactors=F)
  string.inter <- string.all[string.all$experimental>=1 | string.all$database>=1,]
  rm(string.all)

  string.nodes <- unique(c(string.inter[,1], string.inter[,2]))
  string.db <- graphNEL(nodes=string.nodes)
  string.db <- addEdge(string.inter[,1], 
                       string.inter[,2], 
                       string.db)

  expr = read.csv("data/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct", 
                        sep="\t", 
                        skip=2, 
                        stringsAsFactors=F)

  expressed = expr[expr[,"Whole.Blood"] > 0.1, "Name"]
  expressed = sapply(strsplit(expressed, "." ,fixed=T) ,"[", 1)
  library(Homo.sapiens)
  expressed.symbols = select(Homo.sapiens, keys=expressed, keytype="ENSEMBL", columns="SYMBOL")
  expressed.symbols = unique(expressed.symbols[,"SYMBOL"])

  string.nodes = intersect(nodes(string.db), expressed.symbols)
  string.db = subGraph(string.nodes, string.db)

  # get largest connected component
  ig = graph_from_graphnel(string.db)
  cl = clusters(ig)
  keep = nodes(string.db)[cl$membership == which.max(cl$csize)]
  string.db = subGraph(keep, string.db)

  save(file=fcache, string.db)
  } else {
    load(fcache)
  }
  # set to global environment
  assign("STRING.DB", string.db, envir=.GlobalEnv);

  cat("Done.\n");
}

#' Gets ranges in one object close by a set of other ranges
#'
#' Gets the first ranges in subject, which are up-/down-stream and overlapping
#' a range in the query
#'
#' @param query ranges for which to get nearby genes
#' @param subject ranges in which to look for nearby genes
#'
#' @return Either the idx of the hits in the subject if idxs=T, or the identified ranges (idxs=F)
#'
get.nearby.ranges <- function(query, subject) {

  nearby <- function(q,s){
    # get preceding, following and ovberlapping instances of any range in query within subject ranges
    pre <- precede(q, s, select="all")
    fol <- follow(q, s, select="all")
    ove <- findOverlaps(q, s, select="all")
    # combine hits
    h <- unique(c(subjectHits(pre), subjectHits(fol), subjectHits(ove)))

    return(list(hits=h, ranges=s[h]))
  }

  #return a list, where each query gets its nearby ranges annotated with their distance
  res <- lapply(query, function(q) {
    n <- nearby(q, subject);
    n$ranges$distance <- rep(-1, times=length(n$ranges));
    for(i in 1:length(n$ranges)) {
      d <- distance(q,n$ranges[i]);
      n$ranges[i]$distance <- d;
    }
    return(n$ranges);
  });
  return(res);
}

#' Quantile normalization
#'
#' @param x ngenes x nsamples matrix to be normalized
#' @return quantile normalized matrix
#' @export
normalize.quantile <- function(x) {
  x = as.matrix(x)
  o = apply(x, 2, order)
  xsort = x
  for (col in 1:ncol(x)) {
    xsort[,col] = x[o[,col], col]
  }
  means = apply(xsort, 1, mean)
  normalized = matrix(NA, nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
  for (col in 1:ncol(x)) {
    normalized[o[,col], col] = means
  }
  return(normalized)
}


#' Annotate positions with their epigenetic states
#'
#' @param cpg.ranges GRanges object with the positions to annotate
#' @param ids character vector with the roadmap epigenome ids to use
#' @param dir the directory where chromHMM files are stored
#'
#' @return character matrix with length(cpg.ranges) rows and length(ids) columns
#'         with the state annotation of the range in each epigenome
#'
#' @export
chromHMM.annotation <- function(cpg.ranges, 
                                ids, 
                                dir="data/current/roadmap/chromHMM/15state/", 
                                suffix="_15_coreMarks_mnemonics.bed.bgz") {
  require(Rsamtools)
  
  annotation = sapply(ids, function(id) {
    file = paste0(dir, id, suffix)
    avail = as.logical(seqnames(cpg.ranges) %in% seqnamesTabix(file))
    avail.ann = scanTabix(file, param=cpg.ranges[avail])
    avail.ann = sapply(avail.ann, function(x) strsplit(x, "\t")[[1]][4])
    ann = rep(NA, length(cpg.ranges))
    ann[avail] = avail.ann
    return(ann)
  })
  colnames(annotation) = ids
  rownames(annotation) = names(cpg.ranges)

  return(annotation)
}

#' For a set of symbols, gets the emsembl-ids of the corresponding genes from
#' the AnnotationHub. 
#' Note: Currently only implemented for human genes
#'
#' @return A data.frame containng the SYMBOL in one column and the ensemble
#' id in another column
#' 
#' @author Johann Hawe
#' 
#' @date 2017/03/27
#' 
ensembl.from.symbol <- function(symbols, na.drop=T, org="hs"){
  result <- NULL
  if("hs" %in% org){
    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    result <- select(hs, 
                     keytype="SYMBOL", 
                     keys=c(symbols), 
                     columns=c("SYMBOL", "ENSEMBL"))
  } else {
    warning("Organism not supported yet: ", org)
  }
  if(na.drop){
    result <- result[!is.na(result$ENSEMBL),,drop=F]
  }
     
  return(result)
}

#' For a set of ensemble ids, gets the symbols of the corresponding genes from
#' the AnnotationHub. 
#' Note: Currently only implemented for human genes
#'
#' @return A data.frame containng the SYMBOL in one column and the ENSEMBLE
#' id in another column
#' 
#' @author Johann Hawe
#' 
#' @date 2017/03/28
#' 
symbol.from.ensembl <- function(ensemblIDs, na.drop=T, org="hs"){
  result <- NULL
  if("hs" %in% org){
    library(org.Hs.eg.db)
    hs <- org.Hs.eg.db
    result <- select(hs, 
                     keytype="ENSEMBL", 
                     keys=c(ensemblIDs), 
                     columns=c("SYMBOL", "ENSEMBL"))
  } else {
    warning("Organism not supported yet: ", org)
  }
  if(na.drop){
    result <- result[!is.na(result$ENSEMBL),,drop=F]
  }
  # only get the first match for each id
  return(result[!duplicated(result$ENSEMBL),])
}

#' Plot a simple heatmap
#'
#' Uses 'pheatmap' package to plot a very simple heatmap.
#'
#' @param x The data matrix for which to plot the heatmap
#' @param cannot vector of column annotations (must have same names as colnames
#' of x)
#' @param cluster Flag whether to cluster the data in the heatmap or not. default:F
#' @param cors Flag whether correlations are contained in the input matrix x (to
#' appropriately choose the color palette)
#' 
#' @return Grob object returned by pheatmap in the 4th position
#'
simple.heatmap <- function(x, cors=F, cannot=NA, cluster=F) {
  library(pheatmap)
  library(RColorBrewer)
  
  if(cors){
    breaksList = seq(-1, 1, by = 0.1)
  } else {
    breaksList = seq(min(x), max(x), by = 0.1)
  }
  
  cols <- rev(brewer.pal(n = 7, name = "RdYlBu"))
  cols <- colorRampPalette(cols)(length(breaksList))
  
  return(pheatmap(x, annotation_col=cannot, 
           cluster_rows=cluster, 
           cluster_cols=cluster, 
           cex=0.7, color = cols, breaks = breaksList)[[4]])
}


#' For a list of symbols, gets all annotated array ids from the illuminaHumanV3.db
#'
#' @param symbols List of symbols for which to get ids
#' @param mapping Flag whether to return a mapping ID->SYMBOL (default: F)
#' @param as.list Flag whether to return result as list or vector in case mapping=F. (default:F)
#'
#' @author Johann Hawe
#'
probes.from.symbols <- function(symbols, mapping=F, as.list=F, cache.global=F) {
  require(illuminaHumanv3.db)
  if(cache.global & exists("CACHE.ANNOT")){
    annot <- CACHE.ANNOT;
  }
  if(!exists("annot")) {
    annot <- as.list(illuminaHumanv3ALIAS2PROBE)
    annot <- annot[!is.na(annot)]
  }
  if(cache.global & !exists("CACHE.ANNOT")) {
    assign("CACHE.ANNOT", annot, .GlobalEnv)
  }
  annot <- annot[which(names(annot) %in% symbols)];
  if(length(annot) > 0){
    if(mapping) {
      probes <- unlist(annot);
      uprobes <- unique(probes);
      if(length(uprobes) != length(probes)){
        dprobes <- unname(probes[duplicated(probes)])
        warning(paste0("Caution: Some probes had more than one gene annotated,
                       using only one symbol for those probes:\n",
                       dprobes));
      }
      map <- matrix(nrow=length(unlist(annot)))
      rownames(map) <- unlist(annot)
      temp <- lapply(symbols, function(s) {
        ids <- annot[[s]];
        if(!is.null(ids)) {
          map[ids,1] <<- s;
        }
      });
      map <- map[!is.na(map[,1]),,drop=F]
      colnames(map) <- "symbol"
      dropped <- length(symbols) - nrow(map);
      if(dropped>0) {
        cat("Dropped",dropped,"symbols since they had no probe.id available.\n");
      }
      return(map);
    }
    if(as.list) {
      return(annot[symbols])
    } else {
      return(unlist(annot[symbols]))
    }
  } else {
	return(NULL);
  }
}

#' Gets the correlation pvalue matrix between all entities of the supplied
#' data matrix.
#'
#' @param data The data matrix with the individuals in the rows and the measurements
#' in the columns
#' 
#' @param qvalue Flag whether to return qvalues instead of pvalues in the matrix.
#' Default: False
#' 
#' @return A matrix of p-values or q-values between each possible pair
#' of data entities
#'
#' @author Johann Hawe
#' 
get.correlation.pval.matrix <- function(data, qvalue=F){
  library(reshape)
  N <- ncol(data)
  result <- matrix(-1, ncol=N, nrow=N)
  rownames(result) <- colnames(result) <- colnames(data)
  diag(result) <- 1
  for(i in 1:(ncol(data)-1)){
    for(j in (i+1):ncol(data)){
      pval <- cor.test(data[,i], data[,j])$p.value
      result[i,j] <- pval
    }
  }
  if(qvalue){
    # replace pvalues in matrix with the qvalues
    melted <- melt(result)
    # remove self tests
    melted <- melted[which(melted$X1 != melted$X2),]
    melted <- melted[!melted$value==-1,]
    melted$q <- p.adjust(melted$value, "BH")
    for(i in 1:nrow(melted)){
      result[melted[i,"X1"], melted[i,"X2"]] <- melted[i,"q"]
      result[melted[i,"X2"], melted[i,"X1"]] <- melted[i,"q"]
    }
  } else {
    # copy upper tri matrix to lower tri (automatically done when using qvalues)
    for(i in 1:nrow(result)) {
      for(j in 1:(i-1)) {
        result[i,j] <- result[j,i];
      }
    }
  }
  return(result)
}

#' From a symmetric matrix of (correlation) pvalues of measurements, 
#' creates a a graphNEL structure for a given cutoff.
#' 
#' @param cor.pvals The matrix of correlation pvalues
#' @param cutoff The pvalue cutoff to be used
#' 
#' @return A graphNEL object containing significant associations between the 
#' entities
#' 
#' @author Johann Hawe
#' 
create.correlation.network <- function(cor.pvals, cutoff) {
  library(graph)
  
  if(!isSymmetric(cor.pvals)){
    stop("Matrix of correlation pvalues must be symmetric!")
  }
  edges <- matrix(ncol=3, nrow=0)
  colnames(edges) <- c("n1", "n2", "pval")
  
  cn <- colnames(cor.pvals)
  N <- ncol(cor.pvals)
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      pv <- cor.pvals[i,j] 
      if(pv < cutoff){
        edges <- rbind(edges, c(cn[i], cn[j], pv))
      }
    }
  }
  graph <- graphNEL(union(edges[,1], edges[,2]),edgemode = "undirected")
  graph <- addEdge(edges[,1], 
                         edges[,2], 
                         weights=1-as.numeric(edges[,3]),
                         graph)
  return(graph)
}

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
    }
    else {
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

#' Creates the graph at which to start the BDgraph algorithm at
#' 
#' Since it is an MCMC algorithm, starting with this graph will make it more
#' likely to start not too far off of the target distribution/graph.
#' The creation of the start graph is based on the STRING database. 
#' Additionally, predefined edges will be added to the graph.
#'
#' @param nodes All nodes to be incorporated in the graph
#' @param ranges A set of ranges from which cpg and cpg.genes can inferred
#' 
#' @return An incidence matrix of size length(nodes)^2 where 1s indicated edges
#' and 0s indicated absence of edges
#' 
#' @author Johann Hawe
#'
get.g.start <- function(nodes, ranges){
  library(graph)
  
  genes <- nodes[!grepl("^cg|^rs", nodes)]
  
  # create output matrix
  p <- length(nodes)
  out <- matrix(0, nrow=p, ncol=p)
  rownames(out) <- colnames(out) <- nodes
  
  # add STRING connections
  load.string.db()
 
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), genes), STRING.DB)
  sn <- nodes(STRING.SUB)
  em <- t(edgeMatrix(STRING.SUB))
  em <- cbind(sn[em[,1]],sn[em[,2]])
  
  out[em[,1],em[,2]] <- 1
  out[em[,2],em[,1]] <- 1
  
  # create edge matrix for all cpg.gene-cpg edges to be added
  # to g.start
  em <- matrix(ncol=2,nrow=0)
  cpg.genes <- ranges$cpg.genes
  cpg.genes <- cpg.genes[cpg.genes$SYMBOL %in% nodes]
  cpgs <- (ranges$cpgs)
  cpgs <- cpgs[names(cpgs) %in% nodes]
  for(i in 1:length(cpgs)){
    cpg <- cpgs[i]
    g <- get.nearby.ranges(cpgs[i],cpg.genes)[[1]]
    for(row in 1:length(g)){
      em <- rbind(em, 
                  c(names(cpg), 
                    g[row]$SYMBOL))  
    }
    
  }
  
  out[em[,1],em[,2]] <- 1
  out[em[,2],em[,1]] <- 1
  return(out)
}

#' Similar to get.g.start, but simply uses the prior matrix and extracts all
#' pairs with prior>min(priors) to create the g.start
#'
#' @param priors Matrix of prior values for all possible node edges
#' 
#' @author Johann Hawe
#'
get.g.start.from.priors <- function(priors){
  out <- priors  
  out[which(out==min(priors))] <- 0
  out[which(out>min(priors))] <- 1
  return(out)
}

#' Calculate peer factors for given data and covariates
#' 
#' @param data nxg matrix (n=samples, g=genes/variables)
#' @param covariates nxc matrix (n=samples, c=covariates)
#' @param get.residuals Flag whether to directly return the residuals
#' calculated on the data matrix instead of the peer factors calculated
#' @param Nk Number of factors to estimate. Default: N/4
#' 
#' @return Matrix of peer factors calculated
#' 
#' @author Johann Hawe
#' 
#' @date 20170328
#' 
get.peer.factors <- function(data, 
                             covariates=NULL, 
                             get.residuals=F,
                             Nk=ceiling(nrow(data)*0.25)) {
     
  library(peer)
  
  # create model
  model <- PEER();
  # input has to be a matrix!
  PEER_setPhenoMean(model, as.matrix(data));
  
  # add the mean estimation as default since it is recommended in the tutorial
  # of peer. will return Nk+1 factors
  PEER_setAdd_mean(model, TRUE)
  
  # set number of hidden factors to identify. If unknown, a good measure is N/4 (see howto)
  PEER_setNk(model, Nk);
  
  # should not be neccessary but increase anyways
  PEER_setNmax_iterations(model, 5000);
  
  if(!is.null(covariates)){
    # set matrix of known and important covariates, 
    # since we want to acknowledge their effect
    PEER_setCovariates(model, as.matrix(covariates));
  }
  
  # learn 
  PEER_update(model);

  # directly return the residuals if wanted
  if(get.residuals) {
    re <- PEER_getResiduals(model)
    colnames(re) <- colnames(data)
    rownames(re) <- rownames(data)
    return(re)
  }
  
  # get identified factors, contains design.matrix in the first few columns!
  factors <- PEER_getX(model);

  if(!is.null(covariates)){
    # return only the calculated factors, ignoring the original design components  
    factors <- factors[,-c(1:ncol(covariates))]
  }
  
  colnames(factors) <- paste0("f", seq(1:ncol(factors)))
  return(factors);
}

