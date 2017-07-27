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
#' @paran result The bdgraph result
#' @param cutoff Cutoff to be used for posterior probability of edges
#' 
#' @value graph-nel object created from the bdgraph result
#'
#' @author Johann Hawe
#'
graphNEL.from.result <- function(result, cutoff, ranges){
  library(BDgraph)
  library(graph)
  library(igraph)
  
  g.adj <- BDgraph::select(result, cut = cutoff)
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
plot.ggm <- function(graph, id, ranges, dot.out=NULL){
  library(graph)
  library(Rgraphviz)
  library(GenomicRanges)
  
  # remove any unconnected nodes (primarily for bdgraph result, since
  # such nodes are already removed for genenet)
  if(any(graph::degree(graph) == 0)){
    graph <- removeNode(names(which(graph::degree(graph) == 0)), graph)
  }
  
  # add sentinel to network if its is not in there yet or it has been removed
  if(!(id %in% nodes(graph))) {
    graph <- graph::addNode(c(id), graph)
  }
  
  # get trans and cpg gene symbols
  snp.genes <- unique(ranges$snp.genes$SYMBOL)
  cpg.genes <- unique(ranges$cpg.genes$SYMBOL)
  tfs <- NULL
  if("enriched.tfs" %in% names(ranges)){
    tfs <-unique(ranges$enriched.tfs$SYMBOL)
  }
  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=10, 
                          style="filled", fontname="helvetica"), 
                graph=list(overlap="false", root=id, outputorder="edgesfirst"))

  shape = rep("ellipse", numNodes(graph))
  names(shape) = nodes(graph)
  shape[grep("^cg", nodes(graph))] = "box"
  shape[grep("^rs", nodes(graph))] = "box"

  width = rep(0.8, numNodes(graph))
  names(width) = nodes(graph)
  width[grep("cg", nodes(graph))] = 0.4

  height = rep(0.8, numNodes(graph))
  names(height) = nodes(graph)
  height[grep("cg", nodes(graph))] = 0.4

  label = nodes(graph)
  names(label) = nodes(graph)
  label[grep("cg", nodes(graph))] = ""

  col = rep("#ffffff", numNodes(graph))
  names(col) = nodes(graph)
  col[grep("^rs", nodes(graph))] = "#ff0000";
  col[grep("^cg", nodes(graph))] = "#00ff00";
  if(!is.null(tfs)){
    col[tfs] = "green"
  }
  
  penwidth = rep(1, numNodes(graph))
  names(penwidth) = nodes(graph)
  penwidth[snp.genes] = 3
  penwidth[cpg.genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3  
  }
 
  bordercol = rep("black", numNodes(graph));
  names(bordercol) = nodes(graph);
  bordercol[snp.genes] = "red";
  bordercol[cpg.genes] = "green";
 
  bordercol[id] = "red";

  nAttrs = list(shape=shape, label=label, width=width, 
                height=height, penwidth=penwidth, fillcolor=col, 
                color=bordercol)

  ecol = rep("black", numEdges(graph))
  names(ecol) = edgeNames(graph)
  for(i in snp.genes) {
    # color any edge from a SNP to one of our snp genes red
    ecol[grepl("^rs|~rs", names(ecol)) & grepl(i, names(ecol))] = "red"
  }

  dir = rep("none", numEdges(graph))
  names(dir) = edgeNames(graph)

  eAttrs = list(color=ecol, dir=dir)
 
  plot(graph, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)

  if(!is.null(dot.out)){
    # output the dot-information
    #    toDot(graph, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs, subGList=sgs)
    toDot(graph, dot.out, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs)
  }
  
  # return the list of created plotattributes and the possibly modified graph 
  # object
  return(list(graph=graph, nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs))
}

#' Get  STRING interactions (only experimental and db supported)
#'
#' Loads the STRING db saved on disc to a graph R object.
#' Uses a caching mechanism for faster loading of the graph information
#'
#' @return nothing
#'
load.string.db <- function() {
  library(data.table)
  library(graph)
  
  cat("Loading string db.\n")
  fstring <- paste0("results/string.validated.RData")

  cat("Using", fstring, "as cache file.\n")
  
  if(file.exists(fstring)) {
    cat("Loading previously saved network.\n")
    load(fstring)
  } else {
    # load db anew
    string.all <- fread(paste0("data/string/human_gene_hgnc_symbol.links.detailed.v9.0.txt"))
    string.inter <- string.all[experimental>=1 | database>=1]
    rm(string.all)

    STRING.DB <- graphNEL(unique(unlist(string.inter[,1:2, with=F])))
    STRING.DB <- addEdge(unlist(string.inter[,1,with=F]), unlist(string.inter[,2,with=F]), STRING.DB)
    STRING.NODES <- nodes(STRING.DB)
    STRING.EDGES <- edges(STRING.DB)

    save(file=fstring, STRING.DB, STRING.NODES, STRING.EDGES)
  }
  # set vars to global environment
  assign("STRING.DB", STRING.DB, envir=.GlobalEnv);
  assign("STRING.NODES", STRING.NODES, envir=.GlobalEnv);
  assign("STRING.EDGES", STRING.EDGES, envir=.GlobalEnv);

  cat("Done.\n");
}

#' Gets ranges in one object close by a set of other ranges
#'
#' Gets the first ranges in subject, which are up-/down-stream and overlapping
#' a range in the query
#'
#' @param query ranges for which to get nearby genes
#' @param subject ranges in which to look for nearby genes
#' @param idxs flag whether to return the indices of the hits in subject only. default: false
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

#' Gets residuals based on linear models from a data matrix (individuals in the rows)
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
  res <- lapply(colnames(data), function(n) {      
    fm <- as.formula(paste0("`", n, "`~",
                      "1+age+sex"))     
    d <- cbind(data[,n],cov)
    colnames(d)[1] <- n
    return(lm(fm,data=d))
  });

  # build the full residual matrix from model results
  residual.mat <- matrix(data=unlist(lapply(res, resid)), nrow=nrow(data))
  colnames(residual.mat) <- colnames(data);
  rownames(residual.mat) <- rownames(data);

  return(residual.mat)
}

