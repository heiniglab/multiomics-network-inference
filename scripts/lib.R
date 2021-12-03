# ------------------------------------------------------------------------------
#' Method to quickly set a vector of default plotting colors
#' to be used in all generated plots for consistency
#'
#' TODO: parameter for number of colors to get. Then probably
#' the palettes have to be switched accordingly
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------
get_defaultcolors <- function(n=5, name=c("rcb")) {
  library(RColorBrewer)
  # default to rcolorbrewer palette
  cols <- brewer.pal(n, "Set2")
  
  return(cols)
}

# ------------------------------------------------------------------------------
#' Get nodes by type
#'
#' Convenience function to retrieve nodes of a certain type
#' @param g graphNEL locus graph object
#' @param type character indicating the node type
#'
#' @return character vector with the nodes of given type
# ------------------------------------------------------------------------------
get_nodes_by_type <- function(g, type) {
  require(graph)
  n = nodes(g)
  return(n[unlist(nodeData(g, n, type))])
}

# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
#' Get human gene annotation with symbols.
#'
#' This will get the gene annotation annotated with respective gene SYMBOL from UCSC for the hg19 build
#' TODO: allow other builds as well, currently this method is not very flexible...
#'
#' @param drop.nas Flags whether to drop any annotations where the SYMBOL turns out to be NA
#'
#' @return GRanges object, where names(obj) are entrez-ids, and obj$SYMBOLS contains respective gene symbols
#'
#' DEPRECATED
#' Use the new load_gene_annotation() function instead
# ------------------------------------------------------------------------------
get.gene.annotation <- function(drop.nas=TRUE, version=19) {
  warning("DEPRECATED: use file based load_gene_annotation() function instead")
  if(version == 19 | version == 37) {
    library(Homo.sapiens)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    ucsc2symbol = AnnotationDbi::select(Homo.sapiens, keys=keys(Homo.sapiens, keytype="GENEID"),
                                        keytype="GENEID", columns="SYMBOL")
    ga = genes(txdb)
    ga$SYMBOL <- ucsc2symbol[match(values(ga)[,"gene_id"],
                                                ucsc2symbol[,"GENEID"]),"SYMBOL"]
    if(drop.nas){
      ga <- ga[!is.na(ga$SYMBOL)]
    }
  } else if(version == 38) {
    library(annotables)
    grch38 <- data.table(grch38)
#    grch38 <- grch38[grch38$biotype=="protein_coding",]
    ga <- with(grch38,
               GRanges(paste0("chr", chr),
                       IRanges(start, end), strand=strand, SYMBOL=symbol))
  } else {
    stop("Genome annotation version not supported.")
  }

  return(ga)
}

# ------------------------------------------------------------------------------
#' Load gene annotation from a GFF file
#'
#' @param fgene_annot The gene annotation file (GFF format). The method expects
#' to find gene_id, gene_name and gene_biotype in the attributes as well as
#' a single row per gene
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
load_gene_annotation <- function(fgene_annot) {
  require(data.table)
  require(GenomicRanges)

  # load gene annotation
  ga <- fread(fgene_annot)

  # file format is: chr origin type start stop U strand U add_info
  colnames(ga) <- c("chr", "origin", "type", "start", "stop", "score", "strand",
                    "frame", "info")
  # extract ranges
  ra <- with(ga, GRanges(chr, IRanges(start, stop), strand))

  # extract the additional attributes and merge with ranges object
  attrs <- strsplit(ga$info, ";")
  gene_id <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_id",x)], " "), "[[", 2) })
  gene_name <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_name",x)], " "), "[[", 3) })
  gene_biotype <- sapply(attrs, function(x) { sapply(strsplit(x[grepl("gene_type",x)], " "), "[[", 3) })

  # remove any lingering quotes
  gene_id <- gsub("\"", "", gene_id)
  gene_name <- gsub("\"", "", gene_name)
  gene_biotype <- gsub("\"", "", gene_biotype)

  # add to ranges object
  names(ra) <- gene_id
  ra$SYMBOL <- gene_name
  ra$BIOTYPE <- gene_biotype

  # finally, filter out 'misc_RNA' types and unusual chromosomes
  ra <- ra[ra$BIOTYPE != "misc_RNA"]
  ra <- keepStandardChromosomes(ra)

  return(ra)
}

# ------------------------------------------------------------------------------
#'
#' Define a convenience function to get linear model p-values
#'
#' @param modelobject The linear model for which to get the p-value
#'
#' @return The p-value to the given linear model
#'
# ------------------------------------------------------------------------------
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# ------------------------------------------------------------------------------
#' Gets a single snp range from the gtex snp database
#'
#' @author Johann Hawe
#'
# ------------------------------------------------------------------------------
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

#' -----------------------------------------------------------------------------
#' General method to handle either CpG or TSS TFBS context
#'
#' @param entities IDs of CpGs or names of genes for which to get TFBS contet
#' @param fcontext The file containing the precalculated TFBS context. Must
#' match the entity types!
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' -----------------------------------------------------------------------------
get_tfbs_context <- function(entities, fcontext) {
  cont <- readRDS(fcontext)
  tfbs_ann <- cont[rownames(cont) %in% entities,,drop=F]

  if(nrow(tfbs_ann) < 1) return(tfbs_ann)

  entities <- entities[entities %in% rownames(cont)]

  return(tfbs_ann[entities,,drop=F])
}

#' -----------------------------------------------------------------------------
#' Annotates a regulatory graph with appropriate node and
#' edge attributes
#'
#' @param g The graph to be annotated
#' @param ranges The ranges to be used for the annotation. Contains
#' information on which genes are TFs, what the CpG genes are etc.
#' @param ppi_db graphNEL object containing the used PPI database
#' @param fcontext The CpG context file (cpgs annotated with TFBS)
#'
#' @return The same graph instance as given, but annotated with specific
#' node and edge attributes
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
annotate.graph <- function(g, ranges, ppi_db, fcontext){
  # work on all graph nodes
  gn <- nodes(g)

  # check whether we have CpGs
  if(ranges$seed == "meqtl") {
    nodeDataDefaults(g,"cpg") <- F
    nodeData(g, gn, "cpg") <- grepl("^cg", gn)

    nodeDataDefaults(g,"cpg.gene") <- F
    nodeData(g, gn, "cpg.gene") <- gn %in% ranges$cpg_genes$SYMBOL
  } else {
    # eQTLs, we have trans.genes instead of CpGs
    nodeDataDefaults(g,"trans.gene") <- F
    nodeData(g, gn, "trans.gene") <- gn %in% ranges$trans_genes$SYMBOL
  }

  # add common nodeData
  nodeDataDefaults(g,"snp") <- F
  nodeData(g, gn, "snp") <- grepl("^rs", gn)

  nodeDataDefaults(g,"snp.gene") <- F
  nodeData(g, gn, "snp.gene") <- gn %in% ranges$snp_genes$SYMBOL

  nodeDataDefaults(g, "tf") <- F
  nodeData(g, gn, "tf") <- gn %in% ranges$tfs$SYMBOL

  nodeDataDefaults(g, "sp.gene") <- F
  nodeData(g, gn, "sp.gene") <- gn %in% ranges$spath$SYMBOL

  # edge data information
  edgeDataDefaults(g, "isChipSeq") <- FALSE
  edgeDataDefaults(g, "isPPI") <- FALSE

  # add TFBS information
  tfs <- ranges$tfs$SYMBOL
  if(ranges$seed == "meqtl") {
    context <- get_tfbs_context(names(ranges$cpgs), fcontext)
  } else {
    genes <- unique(c(ranges$trans_genes$SYMBOL, ranges$cis_gene$SYMBOL,
                      ranges$tfs$SYMBOL, ranges$snp_genes$SYMBOL, 
                      ranges$spath$SYMBOL))
    context <- get_tfbs_context(genes, fcontext)
  }

  em <- matrix(ncol=2,nrow=0)
  # for all entities
  for(ent in rownames(context)){
    for(tf in tfs) {
      # might be that the TF was measured in more than one cell line
      if(any(context[ent,grepl(tf, colnames(context))]>0)) {
        em <- rbind(em,c(ent,tf))
      }
    }
  }
  em <- filter.edge.matrix(g, em)
  if(nrow(em) > 0){
    edgeData(g, em[,1], em[,2], "isChipSeq") <- T
  }

  # ppi edgedata

  # get subset of edges which are in our current graph
  ss <- subGraph(intersect(nodes(ppi_db), gn), ppi_db)
  edges <- t(edgeMatrix(ss))
  edges <- cbind(gn[edges[,1]], gn[edges[,2]])
  edges <- filter.edge.matrix(g, edges)
  if(nrow(edges) > 0){
    edgeData(g,edges[,1], edges[,2],"isPPI") <- T
  }

  return(g)
}

#' -----------------------------------------------------------------------------
#' Plot a GGM result graph
#'
#' Plots the a built graph (estimated from the sentinel data) using the
#' twopi-visualization. Can be used to retrieve only the plot graph/node/edge
#' attributes when using pdf.out=NULL, dot.out=NULL and plot.on.device=F
#'
#' @param graph: Graph to be plotted (graphNEL)
#'
#' @param id The Id of the sentinel
#' @param dot.out File to which to write the graph in dot format
#'
#' @return List of plot graph attributes and the underlying graph structure
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
plot_ggm <- function(g, id, graph.title=id,
                     plot.on.device=T,
                     dot.out=NULL, ...){

  require(igraph)
  require(graph)
  require(Rgraphviz)

  # get some default colors to be used here
  cols <- get_defaultcolors(n=8)

  # remove any unconnected nodes (primarily for bdgraph result, since
  # such nodes are already removed for genenet)
  if(any(graph::degree(g) == 0)){
    g <- removeNode(names(which(graph::degree(g) == 0)), g)
  }

  # we need an id...
  if(is.null(id)){
    stop("Sentinel ID must not be NULL.")
  }

  # add sentinel to network if it is not in there yet or it has been removed
  if(!is.null(id) && !(id %in% nodes(g))) {
    g <- graph::addNode(c(id), g)
  }

  n <- graph::nodes(g)

  # now get the cluster which contains our sentinel
#  ig = graph_from_graphnel(g)
#  cl = clusters(ig)
#  sentinel_cluster <- cl$membership[names(cl$membership) == id]
#  keep = nodes(g)[cl$membership == sentinel_cluster]
#  g = subGraph(keep, g)
  
  # the remaining nodes
#  n <- keep

  # get trans and cpg gene symbols
  snp.genes <- n[unlist(nodeData(g,n,"snp.gene"))]

  # this is only a quick hack: we plot cpg.genes/trans.genes from the meqtl/eqtl
  # based analyses the same
  if("cpg.gene" %in% names(nodeData(g, nodes(g)[1])[[1]])) {
    assoc_genes <- n[unlist(nodeData(g,n,"cpg.gene"))]
  } else if("trans.gene" %in% names(nodeData(g, nodes(g)[1])[[1]])){
    assoc_genes <- n[unlist(nodeData(g,n,"trans.gene"))]
  }

  tfs <- n[unlist(nodeData(g,n,"tf"))]

  # prepare plot-layout
  attrs <- list(node=list(fixedsize=TRUE, fontsize=14,
                          style="filled", fontname="helvetica"),
                graph=list(overlap="false", root=id[1], outputorder="edgesfirst",
                           label=graph.title, labelloc="top", labeljust="right"))

  shape = rep("ellipse", numNodes(g))
  names(shape) = n
  shape[grep("^cg", n)] = "box"
  shape[grep(id, n)] = "box"

  width = rep(0.8, numNodes(g))
  names(width) = n
  width[grep("cg", n)] = 0.4

  height = rep(0.3, numNodes(g))
  names(height) = n
  height[grep("cg", n)] = 0.4

  label = n
  names(label) = n
  label[grep("cg", n)] = ""

  style <- rep("filled", numNodes(g))
  names(style) <- n

  col = rep("#ffffff", numNodes(g))
  names(col) = n
  col[grep(id, n)] = cols[1]
  col[grep("^cg", n)] = cols[2]
  col[snp.genes] = cols[4]
  col[assoc_genes] = cols[5]
  if(!is.null(tfs)){
    col[tfs] = cols[3]
    tf_cis <- intersect(tfs, snp.genes)
    if(length(tf_cis) > 0) {
      # set two colors since both TF and snp gene
      col[tf_cis] <- cols[6]
    }
  }

  penwidth = rep(1, numNodes(g))
  names(penwidth) = n
  penwidth[snp.genes] = 3
  penwidth[assoc_genes] = 3
  if(!is.null(tfs)){
    penwidth[tfs] = 3
  }

  bordercol = rep("black", numNodes(g));
  names(bordercol) = n;
  bordercol[assoc_genes] = "#e4d7bc";
  bordercol[id] = "#ffe30f";

  nAttrs = list(shape=shape, label=label, style=style,  width=width,
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
  for(cg in assoc_genes){
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

  if(numEdges(g)>60000){
    warning("Skipping plotting on device due to large amount of edges")
  } else if(plot.on.device) {
    plot(g, "twopi", nodeAttrs=nAttrs, edgeAttrs=eAttrs, attrs=attrs, ...)
    # start with empty plot, then add legend
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

    legend("left", legend = c("SNP","SNP gene","TF","assoc gene", "SNP gene/TF"),
           pch=16, pt.cex=3, cex=1.5, bty='n',
           col = c(cols[1], cols[4], cols[3], cols[5], cols[6]), title="graph legend")
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

#' -----------------------------------------------------------------------------
#' Method to quickly filter an edge matrix for only those edges, which are within
#' a specified graph
#'
#' @author Johann Hawe
#'
#' @date 2017/03/13
#'
#' -----------------------------------------------------------------------------
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

#' -----------------------------------------------------------------------------
#' Quantile normalization
#'
#' @param x ngenes x nsamples matrix to be normalized
#' @return quantile normalized matrix
#' @export
#' -----------------------------------------------------------------------------
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
    result <- AnnotationDbi::select(hs,
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
#' @param as_list Flag whether to return result as list or vector in
#' case mapping=F. (default:F)
#'
#' @author Johann Hawe
#'
probes.from.symbols <- function(symbols, mapping=F, as_list=F, annot=NULL) {

  if(is.null(annot)) {
    library(illuminaHumanv3.db)
    annot <- as.list(illuminaHumanv3SYMBOLREANNOTATED)
    annot <- annot[!is.na(annot)]
  }

  annot <- annot[which(annot %in% symbols)]
  annot <- unlist(annot)
  
  if(length(annot) > 0){
    if(mapping) {
      probes <- names(annot)
      uprobes <- unique(probes)
      if(length(uprobes) != length(probes)){
        dprobes <- unname(probes[duplicated(probes)])
        warning(paste0("Caution: Some probes had more than one gene annotated,
                       using only one symbol for those probes:\n",
                       dprobes))
      }
      map <- data.matrix(cbind(names(annot), annot))
      map <- map[!is.na(map[,2]),,drop=F]
      colnames(map) <- c("id", "symbol")
      dropped <- length(symbols) - nrow(map)
      if(dropped>0) {
        cat("Dropped", dropped, "symbols since they had no probe.id available.\n")
      }
      return(map)
    }
    if(as_list) {
      tmp <- sapply(symbols, function(s) { names(annot[annot %in% s])})
      names(tmp) <- symbols
      return(tmp)
    } else {
      return(unique(names(annot[annot %in% symbols])))
    }
  } else {
	return(NULL)
  }
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
  STRING_SUB <- subGraph(intersect(nodes(STRING_DB), genes), STRING_DB)
  sn <- nodes(STRING_SUB)
  em <- t(edgeMatrix(STRING_SUB))
  em <- cbind(sn[em[,1]],sn[em[,2]])

  out[em[,1],em[,2]] <- 1
  out[em[,2],em[,1]] <- 1

  # create edge matrix for all cpg.gene-cpg edges to be added
  # to g.start
  em <- matrix(ncol=2,nrow=0)
  cpg.genes <- ranges$cpg_genes
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
#' pairs with prior>0.5 to create the g.start
#'
#' @param priors Matrix of prior values for all possible node edges
#' @param scaled Do we use a scaled (between 0.5 and 1) version of the priors?
#' Default: FALSE
#' 
#' @author Johann Hawe
#'
get_gstart_from_priors <- function(priors, scaled=FALSE){
  out <- priors
  if(scaled) {
    idx <- out>0.75
  } else {
    idx <- out>0.5
  }
  out[idx] <- 1
  out[!idx] <- 0
  return(out)
}

#' -----------------------------------------------------------------------------
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
#' -----------------------------------------------------------------------------
normalize.expression <- function(data) {
  library(preprocessCore)
  library(sva)

  # quantile normalize
  scaled = normalize.quantiles.robust(t(data))
  rownames(scaled) <- colnames(data)
  colnames(scaled) <- rownames(data)

  # transform the scaled counts to std normal per gene
  stdnorm <- function(x) {
    r = rank(x, ties.method="random")
    qnorm(r / (length(x) + 1))
  }

  # gets p x n matrix
  transformed <- apply(scaled, 1, stdnorm)
  pca <- prcomp(transformed)$x[,1:10]

  # remove first 10 pcs
  corrected <- resid(lm(transformed ~ pca))

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

    locus.graph = graph::addNode(new.nodes, locus.graph)

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
  ## Lov??sz, L. (1993). Random walks on graphs. Combinatorics.

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
  annot <- as.list(illuminaHumanv3SYMBOLREANNOTATED)
  avail <- unique(names(annot))
  ids.avail <- ids[which(ids %in% avail)]
  symbols <- unlist(annot[ids.avail])

  if(mapping) {
    return(data.frame(id=ids.avail,symbol=symbols))
  } else {
    return(symbols)
  }
}

# ------------------------------------------------------------------------------
#' Get the largest connected component in a graphNEL object using igraph
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------
get_largest_cc <- function(g) {
  ig = igraph::graph_from_graphnel(g)
  cl = clusters(ig)
  keep <- cl$membership==which.max(cl$csize)
  keep <- names(cl$membership[keep])
  return(graph::subGraph(keep, g))
}

# ------------------------------------------------------------------------------
#'
#' Method to get gene level estimates from expression probes
#' (=rowMeans on respective probe.ids/gene)
#'
#' @param m The expression matrix with probe ids in columns
#' @param symbols The gene symbols to be used for the summarization.
#' Only probe ids belonging to a gene symbol in this list are used. Rest of
#' probes are dropped.
#'
#' @value Matrix of summarized probe values (per gene probe levels)
#'
#' @author Johann Hawe
#'
#' @date 02/06/2017
#'
#' @export
#'
# ------------------------------------------------------------------------------
summarize <- function(m, symbols){
  # prepare annotation
  library(illuminaHumanv3.db)
  annot <- as.list(illuminaHumanv3SYMBOLREANNOTATED)
  annot <- annot[!is.na(annot)]
  n <- lapply(symbols, function(sym) {
    gid <- probes.from.symbols(sym, annot=annot)
    # for LST1 symbol (and possibly others) there is a probe which should likely
    # be only appointed to the NFKBIL1 gene (according to UCSC genome browser).
    # Might be an error in annotation. Remove probe for LST1 manually for now
    if(sym %in% "LST1"){
      gid <- gid[which(!(gid %in% "ILMN_2046344"))]
    }
    if(is.null(gid)){
      return(NULL)
    }
    if(length(gid)>0){
      if(!all(gid %in% colnames(m))){
        gid <- gid[which(gid %in% colnames(m))]
        if(length(gid)<1){
          warning("None of the probeids for ",
                  sym,
                  " where available in expression data.")
          return(NULL)
        } else {
          warning("Some probes where not available in expression data.")
        }
      }

      g <- m[,gid,drop=F]
      s <- as.matrix(rowMeans(g), ncol=1)
      colnames(s) <- sym
      return(s)
    }
    return(NULL)
  })
  names(n) <- symbols
  n <- n[!sapply(n,is.null)]
  result <- matrix(unlist(n), ncol=length(n), byrow=F)
  colnames(result) <- names(n)
  rownames(result) <- rownames(m)
  return(result)
}

# ------------------------------------------------------------------------------
#' Scans genotype files for SNPs within the provided genomic ranges
#'
#' @param ranges GRanges object containing the ranges which to scan for SNPs
#' @param tabix_file Path to and including file which contains the genotypes,
#' relative to the base directory
#' @param individuals The individual ids for the genotype data in correct order (i.e.
#' as provided in the main directory)
#' @param filter Optional. A list of individuals which should be preselected.
#' Needs to be a subset of `individuals`. Default: NULL
#' 
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
# ------------------------------------------------------------------------------
scan_snps <- function(ranges, tabix_file , individuals, filter=NULL) {
  
  require(Rsamtools)
  require(data.table)
  
  res <- scanTabix(tabix_file, param=ranges)
  
  # scanTabix will also return 'empty' ranges for which no data was found
  res_not_empty <- res[sapply(res, function(x) !length(x)<1)]
  
  print("Scan done, converting Tabix results to data frame.")
  
  if(length(res_not_empty) < 2) {
    collapsed <- paste0(unlist(res_not_empty), "\n")
  } else {
    collapsed <- paste0(unlist(res_not_empty), collapse="\n")
  }
  
  # define columns for filtering
  all_cols <- c("chr", "name", "pos", "orig", "alt", individuals)
  if(!is.null(filter)) {
    filter_cols <- c("chr", "name", "pos", "orig", "alt", filter)
  } else {
    filter_cols <- all_cols
  }
  
  select <- match(filter_cols,
                  all_cols)
  
  # select argument also determines order of columns
  genos <- unique(fread(text=collapsed, sep="\t", header=FALSE, 
                        col.names=filter_cols,
                        stringsAsFactors=F, data.table=FALSE, 
                        select=select))
  
  print("Read data of dimension:")
  print(dim(genos))
  
  # set rownames
  rownames(genos) <- genos$name
  
  # we only return the actual SNP dosages
  return(genos[,6:ncol(genos)])
}

#' -----------------------------------------------------------------------------
#' Get all GWAS traits which were mapped to a specific SNP. Also checks SNPs in
#' LD of the provided SNP (defined via rsquared cutoff)
#'
#' @param snp The snp rsID for which to get GWAS hits
#' @param gwas_table The loaded gwas table. Needs columns 'DISEASE.TRAIT' and
#' SNPS
#' @param get.ld.snps Whether or not to check for LD snps. Default: TRUE
#' @param ld.rsquared Rsquared cutoff to be used in SNiPA to get LD SNPs.
#' Default: 0.95
#' @param collapse Whether to collapse the result traits in a single string.
#' Default: TRUE
#'
#' @author Johann Hawe
#' -----------------------------------------------------------------------------
get_gwas_traits <- function(snp, gwas_table, 
                            get.ld.snps = TRUE,
                            ld.rsquared=0.95, 
                            collapse=TRUE) {
  if(get.ld.snps) {
    # we might miss a hit in which would be in close LD -> get such possible
    # ids via snipa
    snipa_result <- snipa.get.ld.by.snp(snp, rsquare=ld.rsquared)
    aliases <- setdiff(snipa_result$RSALIAS, NA)
    if(length(aliases) > 0) {
      ld_snps <- c(snipa_result$RSID,
                   unlist(strsplit(aliases, ",")),
                   snp)
    } else {
      ld_snps <- c(snp, snipa_result$RSID)
    }
    ld_snps <- unique(ld_snps)
  } else {
    ld_snps <- unique(snp)
  }

  # gather gwas traits for SNPs within this LD block
  gwas_sub <- subset(gwas_table, SNPS %in% ld_snps)
  print(ld_snps)
  if(nrow(gwas_sub) > 0) {
    if(collapse) {
      paste0(unique(gwas_sub$DISEASE.TRAIT), collapse="|")
    } else {
      unique(gwas_sub$DISEASE.TRAIT)
    }
  } else {
    NA
  }
}

#' -----------------------------------------------------------------------------
#' For a list of SNPs, generates a data.frame containing associated GWAS traits
#'
#' Considers proxy SNPs in LD to the supplied SNPs
#'
#' @param snps Vector of SNP rsIDs
#' @param fgwas The file containing the GWAS catalog (ebi)
#'
#' @author Johann Hawe
#'
#' -----------------------------------------------------------------------------
annotate_snps_with_traits <- function(snps, fgwas, drop.nas=TRUE) {
  require(tidyverse)
  source("scripts/snipe.R")

  # load gwas catalog
  gwas <- read_tsv(fgwas) %>% as_tibble(.name_repair="universal")

  # get all traits by snp. Takes a while since it uses SNiPA to get LD SNPs
  traits_by_snp <- lapply(snps, function(s) get_gwas_traits(s, gwas) )

  # prepare result data frame
  annotations <- do.call(rbind.data.frame, traits_by_snp)
  annotations$snp <- snps
  colnames(annotations) <- c("trait", "snp")

  # remove snps without traits?
  if(drop.nas) {
    annotations <- annotations[complete.cases(annotations),]
  }

  annotations <- annotations[,c("snp", "trait")]

  # all done
  annotations
}

# ------------------------------------------------------------------------------
#' Convert a graphNEL object to a table
# ------------------------------------------------------------------------------
graph2table <- function(g) {
  require(dplyr)
  n <- nodes(g)
  em <- t(edgeMatrix(g))
  em <- tibble(n1 = n[em[,1]], n2 = n[em[,2]])
  em
}

# ------------------------------------------------------------------------------
#' Filters a graph for nodes expressed in gtex whole blood
# ------------------------------------------------------------------------------
filter_expression <- function(fgtex, gene_annotation, graph) {
  expr <- fread(fgtex, stringsAsFactors=F)
  expressed <- unlist(expr[`Whole Blood` > 0.1,"Name"])
  expressed.symbols <- ga[intersect(expressed, names(gene_annotation))]$SYMBOL
  
  nodes_to_keep = intersect(nodes(graph), c(expressed.symbols))
  
  return(subGraph(nodes_to_keep, graph))
}
