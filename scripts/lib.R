#' Gets trans-associated CpGs for a sentinel SNP
#'
#' Get all trans cpgs for an individual sentinel SNP using a "cosmo"
#' object.
#'
#' @param sentinel id of the sentinel to analyze
#' @param trans.meQTL the pruned table of trans-associations
#' @param cosmo the cosmo object containing the individual associations
#' @param cpgs the probe ranges of the cpgs to check
#' @param cosmo.idxs flag whether to return the indices in the cosmo object rather than the ones
#' for the cpg list probe ranges.
#'
#' @author Matthias Heinig
#' 
get.trans.cpgs <- function(sentinel, trans.meQTL, cosmo, cpgs=NULL, cosmo.idxs=F) {
  
  pairs = which(trans.meQTL[,"sentinel.snp"] == sentinel)
  
  trans.snp.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.snp"], sep=""), 
                             ranges=IRanges(trans.meQTL[,"interval.start.snp"], 
                                            trans.meQTL[,"interval.end.snp"]))
  trans.cpg.ranges = GRanges(seqnames=paste("chr", trans.meQTL[,"chr.cpg"], sep=""), 
                             ranges=IRanges(trans.meQTL[,"interval.start.cpg"], 
                                            trans.meQTL[,"interval.end.cpg"]))
  
  pair.snps = GRanges(seqnames=paste("chr", cosmo[,"snp.chr"], sep=""), 
                      ranges=IRanges(cosmo[,"snp.pos"], width=1))
  pair.cpgs = GRanges(seqnames=paste("chr", cosmo[,"cpg.chr"], sep=""), 
                      ranges=IRanges(cosmo[,"cpg.pos"], width=2))
  
  pairs = pair.snps %over% trans.snp.ranges[pairs] & 
    pair.cpgs %over% trans.cpg.ranges[pairs]
  
  if(is.null(cpgs) | cosmo.idxs) {
    return(pairs)
  } else {
    is.meQTL = cpgs %over% pair.cpgs[pairs]
    return(is.meQTL)
  }
}

#' Get human gene annotation with symbols.
#'
#' This will get the gene annotation annotated with respective gene SYMBOL from UCSC for the hg19 build
#' TODO: allow other builds as well, currently this method is not very flexible...
#'
#' @param drop.nas Flags whether to drop any annotations where the SYMBOL turns out to be NA
#'
#' @return GRanges object, where names(obj) are entrez-ids, and obj$SYMBOLS contains respective gene symbols
#'
get.gene.annotation <- function(drop.nas=TRUE) {
  library(Homo.sapiens)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  ucsc2symbol = AnnotationDbi::select(Homo.sapiens, keys=keys(Homo.sapiens, keytype="GENEID"), 
                                      keytype="GENEID", columns="SYMBOL")
  ucsc2ens = AnnotationDbi::select(Homo.sapiens, keys=keys(Homo.sapiens, keytype="GENEID"), 
                                      keytype="GENEID", columns="ENSEMBL")
  GENE.ANNOTATION = genes(txdb)
  GENE.ANNOTATION$SYMBOL <- ucsc2symbol[match(values(GENE.ANNOTATION)[,"gene_id"], 
                                              ucsc2symbol[,"GENEID"]),"SYMBOL"]
  GENE.ANNOTATION$ENSEMBL <- ucsc2ens[match(values(GENE.ANNOTATION)[,"gene_id"], 
                                              ucsc2symbol[,"GENEID"]),"ENSEMBL"]
  if(drop.nas){
    GENE.ANNOTATION <- GENE.ANNOTATION[!is.na(GENE.ANNOTATION$SYMBOL)]
  }
  return(GENE.ANNOTATION)
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
#' @param fcontext RData file containing the precomputed chipseq-context
#' 
#' @author Johann Hawe
#' 
get.chipseq.context <- function(cpgs, fcontext){
  load(fcontext)
  tfbs.ann <- tfbs.ann[rownames(tfbs.ann) %in% cpgs,,drop=F]
  
  return(tfbs.ann[cpgs,,drop=F])
}

#' Creates a graphNEL object from a given bdgraph result for a defined cutoff
#' 
#' @param ggm.fit The bdgraph ggm fit
#' @param ranges The ranges of the entities used for the graph fit
#' @param annotate Flag whether to annotate the graph entities with nodeData [T]
#' 
#' @value graph-nel object created from the bdgraph result
#'
#' @author Johann Hawe
#'
graph.from.fit <- function(ggm.fit, ranges, annotate=T){
  library(BDgraph)
  library(graph)
  library(igraph)
  
  # get the graph instance from the ggm fit
  cutoff <- 0.9
  g.adj <- BDgraph::select(ggm.fit, cut = cutoff)
  g <- as_graphnel(graph.adjacency(g.adj, mode="undirected", diag=F))
  
  if(annotate) {
    # set node and edge attributes
    g <- annotate.graph(g, ranges)
  }
  return(g)
}

#' Annotates a regulatory graph with appropriate node and
#' edge attributes
#'
#' @param g The graph to be annotated
#' @param ranges The ranges to be used for the annotation. Contains
#' information on which genes are TFs, what the CpG genes are etc.
#'
#' @return The same graph instance as given, but annotated with specific
#' node and edge attributes 
#' 
#' @author Johann Hawe
#' 
annotate.graph <- function(g, ranges){
  # work on all graph nodes
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
  nodeData(g, gn, "tf") <- gn %in% ranges$tfs$SYMBOL
  
  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% ranges$spath$SYMBOL
  
  # edge data information
  edgeDataDefaults(g, "isChipSeq") <- FALSE
  edgeDataDefaults(g, "isPPI") <- FALSE 
  
  # we will also set tf2cpg priors at this point, so get all TFs
  tfs <- ranges$tfs$SYMBOL
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
  em <- filter.edge.matrix(g,em)
  if(nrow(em) > 0){
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
  if(!is.matrix(em)){
    stop("Provided edge matrix is not an actual matrix.")
  }
  if(nrow(em)<1){
    return(em)
  }
  
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
    pre <- precede(q, s, select="all", ignore.strand=T)
    fol <- follow(q, s, select="all", ignore.strand=T)
    ove <- findOverlaps(q, s, select="all", ignore.strand=T)
    # combine hits
    h <- unique(c(subjectHits(pre), subjectHits(fol), subjectHits(ove)))

    return(list(hits=h, ranges=s[h]))
  }

  #return a list, where each query gets its nearby ranges annotated with their distance
  res <- lapply(query, function(q) {
    n <- nearby(q, subject)
    if(length(n$hits)>0){
      n$ranges$distance <- rep(-1, times=length(n$ranges))
      for(i in 1:length(n$ranges)) {
        d <- distance(q,n$ranges[i])
        n$ranges[i]$distance <- d
      }
      return(n$ranges)
    } else {
      return(NULL)
    }
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


#' For a list of symbols, gets all annotated array ids from the 
#' illuminaHumanV3.db
#'
#' @param symbols List of symbols for which to get ids
#' @param mapping Flag whether to return a mapping ID->SYMBOL (default: F)
#' @param as.list Flag whether to return result as list or vector in 
#' case mapping=F. (default:F)
#'
#' @author Johann Hawe
#'
probes.from.symbols <- function(symbols, mapping=F, as.list=F, annot=NULL) {
  
  if(is.null(annot)) {
    library(illuminaHumanv3.db)
    annot <- as.list(illuminaHumanv3ALIAS2PROBE)
    annot <- annot[!is.na(annot)]
  } 
  
  annot <- annot[which(names(annot) %in% symbols)];
  
  if(length(annot) > 0){
    if(mapping) {
      probes <- unlist(annot)
      uprobes <- unique(probes)
      if(length(uprobes) != length(probes)){
        dprobes <- unname(probes[duplicated(probes)])
        warning(paste0("Caution: Some probes had more than one gene annotated,
                       using only one symbol for those probes:\n",
                       dprobes))
      }
      map <- matrix(nrow=length(unlist(annot)))
      rownames(map) <- unlist(annot)
      temp <- lapply(symbols, function(s) {
        ids <- annot[[s]]
        if(!is.null(ids)) {
          map[ids,1] <<- s
        }
      })
      map <- map[!is.na(map[,1]),,drop=F]
      colnames(map) <- "symbol"
      dropped <- length(symbols) - nrow(map)
      if(dropped>0) {
        cat("Dropped",dropped,"symbols since they had no probe.id available.\n")
      }
      return(map)
    }
    if(as.list) {
      tmp <- annot[symbols]
      names(tmp) <- symbols
      return(tmp)
    } else {
      return(unlist(annot[symbols]))
    }
  } else {
	return(NULL)
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

#' Normalize external (geuvadis, gtex) expression data
#'
#' @param data The expression matrix to normalize, expecting
#' log-transformed RPKM values with individuals in the rows and
#' probes/genes in the columns
#'
#' @return The normalized expression matrix
#'
#' @author Johann Hawe
#'
normalize.expression <- function(data) {
  library(preprocessCore)
  library(peer)

  # quantile normalize 
  scaled = normalize.quantiles(t(data))
  rownames(scaled) <- colnames(data)
  colnames(scaled) <- rownames(data)

  # transform the scaled counts to std normal per gene
  stdnorm <- function(x) {
    r = rank(x, ties.method="random")
    qnorm(r / (length(x) + 1))
  }
  transformed <- apply(scaled, 1, stdnorm)

  # remove peer factors
  corrected <- correct.peer(transformed, Nk = 10)
  colnames(corrected) <- colnames(transformed)
  return(corrected)
}

#' Calculate peer factors for given data and covariates
#' Then get residual data matrix 
#'
#' @param data nxg matrix (n=samples, g=genes/variables)
#' @param covariates nxc matrix (n=samples, c=covariates)
#' @param get.residuals Flag whether to directly return the residuals
#' calculated on the data matrix instead of the peer factors calculated
#' @param Nk Number of factors to estimate. Default: N/4
#' 
#' @return Matrix of corrected expression data
#' 
#' @author Johann Hawe
#' 
#' @date 20170328
#' 
correct.peer <- function(data,
                             covariates=NULL,
                             Nk=ceiling(nrow(data)*0.25)) {
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

  re <- PEER_getResiduals(model)
  colnames(re) <- colnames(data)
  rownames(re) <- rownames(data)
  return(re)
}

#' Gets residuals based on linear models from a data matrix (individuals in the rows). This
#' is specifically used for the Gtex expression matrix rather than the cohort data
#'
#' Calculates the residual matrix for a given expression matrix 
#' considering available covariates. Uses linear model. Expects expression probe ids to start
#' with "ILMN_". Supplied covariates must not start with either of the prefixes
#'
#' @param data the matrix for which to calculate the covariates. Needs to contain covariates themselfes plus
#' either the methylation or expression probes for which to get the residuals
#'
#' @return A matrix  where in the colums are the measured entities (e.g.probes) 
#' and in the rows are the samples, containing the calculated residuals. Covariates supplied in the
#' input matrix are discared
#'
get.residuals <- function(data, cov) {
  # process GTEX data
  res <- lapply(colnames(data), function(n) {      
    fm <- as.formula(paste0("`", n, "`~",
                      "1+age+sex+batch1+batch2"))     
    d <- cbind(data[,n],cov)
    colnames(d)[1] <- n
    return(lm(fm,data=d))
  });
  residual.mat <- matrix(data=unlist(lapply(res, resid)), nrow=nrow(data))
  colnames(residual.mat) <- colnames(data);
  rownames(residual.mat) <- rownames(data);
  return(residual.mat)
}

#' Small helper for a list() function automatically setting
#' the variable names as the names of the created list
#' 
#' compare: https://stackoverflow.com/a/21059868
#'
#' @author Johann Hawe
#' 
listN <- function(...){
  ln <- list(...)
  names(ln) <- as.character(substitute(list(...)))[-1]
  ln
}

#' Method to quickly rewire a given graph
#' 
#' @param g The graphNEL object to be rewired
#' @param p The probability with which an edge gets rewired
#' @param ei Dataframe of class edge_info
#' @return The rewired graph object
#' 
rewire_graph <- function(g, p, ei) {
  if(!"edge_info" %in% class(ei)){
    stop("Given data.frame is not an edge info object")
  }
  # do nothing...
  if(p == 0) {
    return(list(gnew=g, to_switch=NULL))
  }
  
  # get node degrees
  degs <- cbind.data.frame(sort(graph::degree(g), decreasing = T))
  colnames(degs) <- c("degree")
  degs$node <- rownames(degs)
  
  # number of edges to rewire
  N <- floor(p * numEdges(g))
  
  # get the 0-degree nodes which do not have any prior associated
  ei_priors <- ei[ei$keep_full,,drop=F]
  nodes <- nodes(g)
  nodes_no_priors <- nodes[!(nodes %in% ei_priors$n1) &
                             !(nodes %in% ei_priors$n2)]
  
  # select nodes to switch such that sum of node degree == N
  idxs <- which(degs$degree > 0)
  # keep within this number of nodes if possible
  Nmax <- length(nodes_no_priors)
  lo <- N-1
  up <- N+1
  best <- NULL
  best_diff <- 0.5
  
  while(best_diff>0) {
    res <- degreesum_in_range(degs, sample(idxs), 
                                       lo, up, 
                                       g)
    if(!is.null(res) & !is.null(res$nodes)) {
      s <- sum(res$nodes$degree)
      a <- res$adjustment
      # get 'best' result -> closest to desired number of change
      d <- abs(s-a-N)
      n <- nrow(res$nodes)
      if(d < best_diff & n <= Nmax) {
        best <- res$nodes
        best_diff <- d
      }
    }
  }
  Ns <- nrow(best)
 
  if(length(nodes_no_priors) > 0) {
    # define switchings
    Nnp <- length(nodes_no_priors)
    
    to_switch <- cbind.data.frame(n_prior=best$node, 
                                  n_prior_degree=best$degree,
                                  n_noprior=sample(nodes_no_priors, 
                                                   size = Ns, 
                                                   replace = Ns>Nnp),
                                  stringsAsFactors=F)
    
    gnew <- switch_nodes(g, to_switch)
   
  } else {
    stop("Not implemented yet: There where no nodes without priors")
  }
  
  return(listN(gnew, to_switch))
}

#' Creates a graph with switched nodes as compared to a base graph
#' 
#' @param g_old
#' @param to_switch
#'
#' @author Johann Hawe
#' 
switch_nodes <- function(g_old, to_switch) {
  # get the full edge matrix and switch nodes
  n <- nodes(g_old)
  em <- t(edgeMatrix(g_old))
  em <- cbind.data.frame(n1=n[em[,1]], 
                         n2=n[em[,2]],
                         stringsAsFactors=F)
  
  for(i in 1:nrow(to_switch)) {
    prior_node <- to_switch[i,"n_prior"]
    noprior_node <- to_switch[i,"n_noprior"]
    # replace all prior node occurrences with the nonprior node
    em[em$n1 == prior_node,"n1"] <- noprior_node
    em[em$n2 == prior_node,"n2"] <- noprior_node
  }
  
  g_new <- graphNEL(n)
  g_new <- addEdge(em[,1], em[,2], g_new)
  return(g_new)
}

#' Gets a set of nodes from a dataframe with annotated degree,
#' where the sum of all degrees is within a certain window.
#' 
#' @param df The dataframe containing the nodes and their degrees (degree column)
#' @param idxs Indicies of the nodes in the df to be used
#' @param lo Lower bound for the sum
#' @param up Upper bound for the sum
#' @param g The graph from which the nodes originate. If given, it is ensured
#' that nodes in the final collection of nodes are not neighbours within the graph
#' 
#' @return Set of nodes for which the sum of degrees is between lo and up
#' 
#' @author Johann Hawe
#' 
degreesum_in_range <- function(df, idxs, lo, up, g) {
  
  # get the  center value for exact matching
  total <- (lo+up)/2
  
  # check special cases
  if(total == 0) {
    return(NULL)
  } else if(total == 100) {
    return(df[idxs,,drop=F])
  }
  
  rand <- c()
  # counter for total running sum of degrees we added
  running_sum <- 0
  # counter of substracted information, i.e. whenever we find a new node has
  # neighbours within our collection, we reduce the total amount of degrees found
  # by the respective number of neighbours. We need this information for 
  # the evaluation of the found set of nodes
  total_subst <- 0
  
  for(i in idxs) {
    pi <- df[i,,drop=F]
    pi <- cbind(pi, idx=i)
    v <- pi[,"degree"]
    s <- running_sum + v
    subst <- 0
    # neighbour already in collection? handle that
    # adjust total sum of degrees we switch when using this one
    if(!is.null(nrow(rand))) {
        nn <- rand[,"node"]
        # get all neighbours of current node to check how many are already in our
        # collection
        cn <- pi[,"node"]
        neigh <- unlist(graph::adj(g, cn))
        subst <- sum(nn %in% neigh)
        s <- s - subst
    }
    # did not reach lower bound yet
    if(s<lo) {
      running_sum <- s
      total_subst <- total_subst + subst
      rand <- rbind(rand,pi)
      next
    }
    # above lower bound
    if(s>=lo) {
      # match exactly? then done
      if(s==total) {
        running_sum <- s
        total_subst <- total_subst + subst
        rand <- rbind(rand,pi)
        break
      }
      # below upper bound
      if(s<=up) {
        rand <- rbind(rand,pi)
        running_sum <- s
        # first above exact? then done
        if(s>total) {
          running_sum <- s
          total_subst <- total_subst + subst
          break
        }
      } else {
        # exceeded with current value, stop
        running_sum <- s
        total_subst <- total_subst + subst
        break
      }
    }
  }
  
  return(list(nodes=rand, adjustment=total_subst))
}

#' Helper method calculating the percentage of edges in the given graph object
#' which have a prior according to the given prior matrix
#' 
#' @param g The graph for which to check the edges
#' @param priors The symmetric matrix of edge priors. The lowest value
#' of this matrix will be considered as the pseudo prior. Only values 
#' larger than the pseudo prior will be considered as 'priors'.
#' 
#' @return Percentage of edges in the graph which have a prior
get_percent_prioredges <- function(g, priors) {
  library(graph)
  
  # get pseudo prior value
  pp <- min(priors)
  
  # get the edgematrix with nodes
  n <- graph::nodes(g)
  
  # sanity check
  if(!all(n %in% colnames(priors))) {
    warning("Not all nodes are in prior matrix.")
    return(NULL)
  }
  em <- t(edgeMatrix(g))
  if(nrow(em) > 0) {
    count <- 0
    em <- cbind.data.frame(n[em[,1]], 
                n[em[,2]],
                stringsAsFactors=F)
    for(i in 1:nrow(em)) {
      v <- priors[em[i,1], em[i,2]]
      # check if value is larger than pseudo prior
      if(v>pp) {
        count <- count + 1
      }
    }
    return(count/nrow(em))
  } else {
    warning("No edges in graph.")
    return(NULL)
  }
}

remove_unconnected_nodes <- function(g) {
  to_remove <- nodes[graph::degree(g,nodes(g))<1]
  if(length(to_remove) > 0) {
    g <- removeNode(to_remove, g)
  }
}

#' Samples a graph from the given prior matrix
#' Uses the colnames of the (symmetric) given prior
#' matrix to determine graph nodes
#' 
#' @param priors The prior matrix
#'
#' @author Johann Hawe
#' 
sample_prior_graph <- function(priors, sentinel) {
  
  library(reshape2)
  
  # assume pseudo.prior as min-value (should always be
  # the case)
  pseudo.prior <- min(priors)
  
  # the nodes we want to operate on
  nodes <- colnames(priors)
  
  # create base graph
  g <- graphNEL(nodes, edgemode = "undirected")
  
  # iterate over each combination of nodes
  temp <- priors
  temp[upper.tri(temp, T)] <- NA
  ee <- melt(temp, na.rm=T,stringsAsFactors=F)
  colnames(ee) <- c("n1", "n2", "prior")
  ee$n1 <- as.character(ee$n1)
  ee$n2 <- as.character(ee$n2)
  
  # create names
  rownames(ee) <- paste(ee$n1, ee$n2, sep="_")
  
  # flag for sampling
  ee$keep <- F
  
  # iterate over each pair and determine whether to add
  # the pair as an 'true' edge
  for(i in 1:nrow(ee)) {
    # get prior
    p <- ee[i,"prior"]
    if(p>pseudo.prior) {
      v <- runif(1)
      if(v<=p) {
        ee[i,"keep"] <- T
      }
    }
  }
  
  # create a full graph for comparison
  ee$keep_full <- ee$prior > pseudo.prior
  gfull <- addEdge(ee[ee$keep_full,1], ee[ee$keep_full,2], g)
  
  # keep only relevent edges and add to graph
  to_add <- ee[ee$keep,,drop=F]
  if(nrow(to_add)>0) {
    # make sure we add the SNP
    # -> if its not in the list of edges to keep,
    # we simply switch one random connection with it
    # if it doesn't have a prior, though, we just ignore it
    if(!(sentinel %in% c(to_add[,"n1"], to_add[,"n2"]))) {
      temp <- priors[,sentinel]
      temp <- temp[temp>pseudo.prior]
      if(!length(temp) == 0) {
        # select random name
        temp2 <- sample(names(temp),1)
        to_switch <- sample(1:nrow(to_add),1)
        old <- to_add[to_switch,]
        old["keep"] <- F
        
        to_add[to_switch,] <- c(sentinel, temp2, priors[sentinel,temp2], T, T)
        rownames(to_add)[to_switch] <- paste0(sort(c(sentinel, temp2)), 
                                              collapse="_")
        ee[ee$n1 == old$n1 & ee$n2 == old$n2,] <- old
      }
    }
  } else {
    # none of the edges made it
    # for now add at least the edges which got any prior so that we do not 
    # get an empty graph here; add random SNP edge to have the SNP available
    warning("Creating full graph as observed graph.")
    ee[ee$prior>pseudo.prior,"keep"] <- T
    ee[which(grepl("^rs", ee$n1) | grepl("^rs", ee$n2))[1], "keep"] <- T
    to_add <- ee[ee$keep,,drop=F]
  }
  
  g <- addEdge(to_add[,1], to_add[,2], g)
  class(ee) <- c(class(ee), "edge_info")
  return(list(sample_graph=g, full_graph=gfull, 
              edge_info=ee, sentinel=sentinel))
}

#' Augment grahs
#'
#' Add TF to CpG edges and trans gene to SNP edges to graphs. This function
#' assumes that you have an object tfbs.ann in the environment. It is not
#' passed as argument (bad style) becuase it is too big (passing by value).
#'
#' @param graphs list of graph objects to add nodes and edges to
#' @param sentinel character id of the sentinel SNP
#' @param trans.genes character vector of genes in the trans locus
#' @param trans.cpgs character vector of CpGs with trans meQTLs
#' @param tfbs.ann logical indicator matrix with CpGs in the rows and
#'        transcription factors as columns
#'
#' @return list of graph objects with the nodes and edges added
#'
#' @author Matthias Heinig
#' 
#' @imports reshape
#' @imports graph
#' @export
add.to.graphs <- function(graphs, sentinel, trans.genes, trans.cpgs, tfbs.ann) {
  require(reshape)
  require(graph)
  
  tf.edges = melt(tfbs.ann[trans.cpgs,])
  colnames(tf.edges) = c("cpg", "condition", "adjacent")
  tf.edges = tf.edges[tf.edges[,"adjacent"],]
  for (i in 1:2) {
    tf.edges[,i] = as.character(tf.edges[,i])
  }
  tf = as.character(sapply(strsplit(tf.edges[,"condition"], ".", fixed=T), 
                           "[", 1))
  tf.edges = data.frame(tf.edges, tf, stringsAsFactors=F)
  tf.edges = tf.edges[!duplicated(paste(tf.edges[,"tf"], tf.edges[,"cpg"])),]
  
  ## put together the graphs
  out = list()
  
  for (graph.idx in 1:length(graphs)) {
    locus.graph = graphs[[graph.idx]]
    
    ## filter for TFs that are in the graph already
    use.tf.edges = tf.edges[tf.edges[,"tf"] %in% nodes(locus.graph),]
    
    new.nodes = unique(c(use.tf.edges[,"tf"], trans.cpgs, 
                         sentinel, trans.genes))
    new.nodes = setdiff(new.nodes, nodes(locus.graph))
    
    locus.graph = addNode(new.nodes, locus.graph)
    
    ## also add some meta data (the type of nodes)
    nodeDataDefaults(locus.graph, "tf") = FALSE
    nodeData(locus.graph, unique(use.tf.edges[,"tf"]), "tf") = TRUE
    
    nodeDataDefaults(locus.graph, "cpg") = FALSE
    nodeData(locus.graph, trans.cpgs, "cpg") = TRUE
    
    nodeDataDefaults(locus.graph, "snp") = FALSE
    nodeData(locus.graph, sentinel, "snp") = TRUE
    
    nodeDataDefaults(locus.graph, "trans.gene") = FALSE
    nodeData(locus.graph, trans.genes, "trans.gene") = TRUE
    
    ## add edges for the tfbs
    locus.graph = addEdge(use.tf.edges[,"tf"], 
                          use.tf.edges[,"cpg"], locus.graph)
    
    ## add edges for the connection of the locus and its genes
    locus.graph = addEdge(rep(sentinel, length(trans.genes)), 
                          trans.genes, locus.graph)
    
    out[[graph.idx]] = locus.graph
  }
  
  return(out)
}

## this function is from destiny and uses the arpack function
eig.decomp <- function(M, n.eigs, sym) {
  n = nrow(M)
  f <- function(x, extra=NULL) as.matrix(M %*% x)
  wh <- if (sym) 'LA' else 'LM'
  ## constraints: n >= ncv > nev
  ar <- arpack(f, sym = sym, options = list(
    which = wh, n = n, ncv = min(n, 4*n.eigs), nev = n.eigs + 1))
  if (!sym) {
    ar$vectors <- Re(ar$vectors)
    ar$values  <- Re(ar$values)
  }
  
  ## check for negative values and flip signs
  neg = which(ar$values < 0)
  for (n in neg) {
    ar$values[n] = -ar$values[n]
    ar$vectors[,n] = -ar$vectors[,n]
  }
  
  return(ar)
}

graph2sparseMatrix <- function(g) {
  em = edgeMatrix(g)
  Asparse = sparseMatrix(em[1,], em[2,], x=1, dims=rep(numNodes(g), 2), 
                         dimnames=list(nodes(g), nodes(g)))
  
  ## make symmetric
  Asparse = Asparse + t(Asparse)
  Asparse[Asparse > 1] = 1
  
  return(Asparse)
}


## approximate the infinite sum
## in the paper and notes this is matrix is called M
## n.eigs gives the number of eigenvectors used for the approximation
## if this is smaller than 2, the pseudo inverse will be computed
## from and to are the nodes for which we need the transition probabilties
## this avoids to compute the full propagation matrix which is to large
## when from and to are NULL the full matrix will be computed
## sum can be "none", "from", "to" or "both" and will sum up the propagation
## matrix over one or the other or both lists
## "none" returns a "from" by "to" matrix
## "from" returns a vector of length nrow(Asparse)
## "to" returns a vector of length nrow(Asparse)
## "both" returns a nrow(Asparse) by 2 matrix with the from and to vectors
propagation <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {
  require(igraph)
  require(Matrix)
  ## transition matrix
  ## transition = Asparse / rowSums(Asparse)
  
  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(Matrix::rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root
  
  if (is.null(from)) {
    from = 1:nrow(Asparse)
  }
  if (is.null(to)) {
    to = 1:nrow(Asparse)
  }
  if (is.logical(from)) {
    from = which(from)
  }
  if (is.logical(to)) {
    to = which(to)
  }
  dnames = list(from, to)
  ## if we have character we need to match
  if (is.character(from)) {
    from = match(from, rownames(Asparse))
  }
  if (is.character(to)) {
    to = match(to, colnames(Asparse))
  }
  
  ## setup different ways of summarizing the propagation matrix
  if (sum == "from") {
    transform = rowSums
    propagation = double(0, length(to))
    names(propagation) = dnames[[2]]
  } else if (sum == "to") {
    transform = colSums
    propagation = double(0, length(from))
    names(propagation) = dnames[[1]]
  } else if (sum == "both") {
    propagation = matrix(0, nrow=nrow(Asparse), ncol=2, 
                         dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    propagation = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }
  
  if (n.eigs > 0) {
    ## approximate using the eigen vectors
    
    if (n.eigs < nrow(Psym)) {
      ## just use a subset of eigen vectors
      eig.M <- eig.decomp(Psym, n.eigs, TRUE)
      
    } else {
      ## use all eigen vectors (only for small matrices)
      eig.M <- eigen(Psym)
    }
    
    ## for both computations we discard the first eigenvector as it
    ## represents the stationary distribution
    if (sum %in% c("none", "from", "to")) {
      for (i in 2:n.eigs) {
        increment = (eig.M$values[i] / (1 - eig.M$values[i])) *
          eig.M$vectors[from,i] %*% t(eig.M$vectors[to,i])
        propagation = propagation + transform(increment)
      }
      
    } else {
      ## when sum is "both" we iterate over the "from" and "to" nodes and add
      ## the from and to vectors incrementally to save memory
      for (i in 2:n.eigs) {
        weight = (eig.M$values[i] / (1 - eig.M$values[i]))
        ## increment the "from" sum vector
        for (ff in from) {
          propagation[,1] = propagation[,1] + weight * (eig.M$vectors[ff,i] * eig.M$vectors[,i])
        }
        ## increment the "to" sum vector
        for (tt in to) {
          propagation[,2] = propagation[,2] + weight * (eig.M$vectors[tt,i] * eig.M$vectors[,i])
        }
      }
    }
  } else {
    ## compute using the pseudo inverse
    
    ## the first eigenvector corresponds to stationary distribution defined
    ## by the degrees of the nodes
    ## Attention: this is actually the first eigenvector of the asymmetric
    ## Psym matrix
    d = rowSums(Asparse)
    v = sum(d)
    phi0 = d / v
    
    ## there is an equivalence of the eigenvectors of the symmetric and
    ## asymetric matrix:
    ## EV(sym) = D^-0.5 EV(asym)
    phi0 =  phi0 / sqrt(d)
    
    ## in the dtp package there is a length normalization step that matches
    ## the vectors exactly (normalizing to unit length vectors)
    phi0 = phi0 / sqrt(sum(phi0^2))
    
    n <- nrow(Psym)
    inv <- solve(Diagonal(n) - Psym + phi0 %*% t(phi0))
    propagation = inv - Diagonal(n)
    propagation = propagation[from, to]
    if (sum %in% c("none", "from", "to")) {
      propagation = transform(propagation)
    } else {
      stop("sum = 'both' not implemented yet for pseudo inverse")
    }
  }
  return(propagation)
}





## alternatively we can compute expected hitting times for the random walks
## n.eigs gives the number of eigenvectors used for the approximation
## if this is smaller than 2, the pseudo inverse will be computed
## from and to are the nodes for which we need the transition probabilties
## this avoids to compute the full propagation matrix which is to large
## when from and to are NULL the full matrix will be computed
## sum can be "none", "from", "to" or "both" and will sum up the propagation
## matrix over one or the other or both lists
## "none" returns a "from" by "to" matrix
## "from" returns a vector of length nrow(Asparse)
## "to" returns a vector of length nrow(Asparse)
## "both" returns a nrow(Asparse) by 2 matrix with the from and to vectors
## the i-th component of the "from sum" represents the average hitting time
## going from node i to any of the "to" nodes
## the i-th component of the "to sum" represents the average hitting time
## going from any "from" nodes to node i
hitting.time <- function(Asparse, n.eigs=20, from=NULL, to=NULL, sum="none") {
  
  ## we need the number of edges and the number of nodes and the degree of nodes
  numEdges = sum(Asparse) / 2 + sum(diag(Asparse))
  numNodes = nrow(Asparse)
  deg = rowSums(Asparse)
  
  ## symmetric transition matrix
  D.root = Diagonal(nrow(Asparse), 1 / sqrt(rowSums(Asparse)))
  Psym = D.root %*% Asparse %*% D.root
  
  if (is.null(from)) {
    from = 1:numNodes
  }
  if (is.null(to)) {
    to = 1:numNodes
  }
  if (is.logical(from)) {
    from = which(from)
  }
  if (is.logical(to)) {
    to = which(to)
  }
  dnames = list(from, to)
  ## if we have character we need to match
  if (is.character(from)) {
    from = match(from, rownames(Asparse))
  }
  if (is.character(to)) {
    to = match(to, colnames(Asparse))
  }
  
  ## setup different ways of summarizing the hitting.time matrix
  if (sum == "from") {
    transform = rowMeans
    hitting.time = double(0, length(to))
    names(hitting.time) = dnames[[2]]
  } else if (sum == "to") {
    transform = colMeans
    hitting.time = double(0, length(from))
    names(hitting.time) = dnames[[1]]
  } else if (sum == "both") {
    hitting.time = matrix(0, nrow=numNodes, ncol=2, dimnames=list(rownames(Asparse), c("from", "to")))
  } else if (sum == "none") {
    transform = function(x) {return(x)}
    hitting.time = matrix(0, nrow=length(from), ncol=length(to), dimnames=dnames)
  }
  
  ## approximate using the eigen vectors
  
  if (n.eigs < nrow(Psym)) {
    ## just use a subset of eigen vectors
    eig.M <- eig.decomp(Psym, n.eigs, TRUE)
    
  } else {
    ## use all eigen vectors (only for small matrices)
    eig.M <- eigen(Psym)
  }
  
  ## Hitting time computed according to Theorem 3.1 from
  ## Lovász, L. (1993). Random walks on graphs. Combinatorics.
  
  compute.hitting.time <- function(from.nodes, to.nodes) {
    sapply(to.nodes, function(tt)
      sapply(from.nodes, function (ff) {
        ht = 2 * numEdges * sum(sapply(2:n.eigs, function(i) {
          v = eig.M$vectors[,i]
          return(1 / (1 - eig.M$values[i]) *
                   (v[tt]^2 / deg[tt] - v[ff] *
                      v[tt] / sqrt(deg[ff] * deg[tt])))
        }))
      }))
  }
  
  if (sum %in% c("none", "from", "to")) {
    for (i in 2:n.eigs) {
      increment = compute.hitting.time(from, to)
      hitting.time = hitting.time + transform(increment)
    }
    
  } else {
    ## increment the "from" sum vector
    for (tt in to) {
      hitting.time[,1] = hitting.time[,1] + compute.hitting.time(1:numNodes, tt)
    }
    hitting.time[,1] = hitting.time[,1] / length(to)
    
    ## increment the "to" sum vector
    for (ff in from) {
      hitting.time[,2] = hitting.time[,2] + compute.hitting.time(ff, 1:numNodes)
    }
    hitting.time[,2] = hitting.time[,2] / length(from)
  }
  
  return(hitting.time)
}


## find the shortest path with minimal node weight
min.node.weight.path <- function(g, weights, from, to) {
  ## the problem can be transformed into a directed graph problem where all
  ## incoming edges are assigned the node weight
  
  require(graph)
  require(RBGL)
  
  h = graphNEL(nodes(g), edgemode="directed")
  em = edgeMatrix(g)
  h = addEdge(nodes(g)[em[1,]], nodes(g)[em[2,]], h, weights[em[2,]])
  h = addEdge(nodes(g)[em[2,]], nodes(g)[em[1,]], h, weights[em[1,]])
  
  return(sp.between(h, from, to))
}

#' Get gene symbols for the given (array) probe ids
#' 
#' @param ids List of probe ids for which to get the symbols
#' @param mapping Whether to return a mapping (probe->symbol); default: F
#' 
#' @author Johann Hawe
#' 
#' @return If mapping=F: All symbols referenced by the given probe ids. Otherwise
#' a mapping of all probeids to the corresponding symbols (a data.frame)
#'
symbols.from.probeids <- function(ids, mapping=F){
  library(illuminaHumanv3.db)
  annot <- illuminaHumanv3SYMBOLREANNOTATED
  avail <- mappedkeys(annot)
  ids.avail <- ids[which(ids %in% avail)]
  
  symbols <- unlist(as.list(annot[ids.avail]))
  
  if(mapping) {
    return(data.frame(id=ids.avail,symbol=symbols))
  } else {
    return(symbols)
  }
}
