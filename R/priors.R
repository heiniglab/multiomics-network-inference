#' Get priors for links between entities
#'
#' Given appropriate information, creates a matrix for all GGM entitites (i.e.
#' snps, cpgs and cis-/trans-genes) reflecting respective prior probabilities for links
#' between them (currently only based on distances of snps/cpgs <-> genes)
#'
#' @param ranges Set of ranges with the annotated data (e.g. snp.ranges, cpgs etc.) for which to
#' extract the link priors
#'
#' @param nodes All the nodes/entities being analyzed.
#' 
#' @return A square matrix of link-priors
#'
get.link.priors <- function(ranges, nodes) {
  
  # assume single SNP
  message("WARNING: Assuming single SNP in givne data.")
  id <- nodes[grepl("^rs", nodes)]
  
  # load string db
  load.string.db();
  # load predefined prior definitions from gtex
  load.gtex.priors(id);
  
  # default window (needed for distance calculation) and number of bins
  wnd <- 1e6;
  nbins <- 200;
   # create bin vector
  step <- wnd/nbins;
  bins <- seq(0,wnd-step, by=step);
  
  
  # get distances to all genes for cpgs and snps
  cpg.dists <- get.nearby.ranges(ranges$cpgs, 
                                 promoters(ranges$cpg.genes))
  all.dists <- append(cpg.dists, snp.dists);

  # sanity check
  if(!all(names(all.dists) %in% nodes)){
    cat("Some elements were not available in data:\n")
    na <- names(all.dists)
    cat(na[which(!na %in% nodes)])
    all.dists <- all.dists[!names(all.dists) %in% na]
    warning("Some data mysteriously went missing..")
  }
  
  # now build the prior matrix using a pseudo prior
  pseudo.prior <- 0.0000001
  priors <- matrix(data = pseudo.prior, nrow = length(nodes), ncol=length(nodes))
  colnames(priors) <- rownames(priors) <- nodes

  # load annotation needed for cpg2gene priors
  epigen.states <- read.table("data/current/epigenetic_state_annotation_weighted_all_sentinels.txt", 
                        header = T, sep="\t", row.names = 1)
  
  ## SET CPG-CPGGENE PRIOR
  for(i in names(all.dists)){
    gene.ranges <- all.dists[[i]];
    temp <- lapply(gene.ranges, function(r) {
      s <- r$SYMBOL;
      if(s %in% colnames(priors)) {
        # set the basic prior based on distance (the larger the distance, the lower the prior should be)
        # but scale it to be between 0 and 1
        dist <- r$distance
        if(!is.na(dist)){
          # if the cpg is in range of the tss (within 200bp), 
          # set a specific prior for active TSS... 
          p <- pseudo.prior
          if(dist <= 200){
            # set the cpg.state prior
            p <- p + sum(epigen.states[i, c("Active.TSS", 
                                        "Flanking.Active.TSS", 
                                        "Bivalent.Poised.TSS",
                                        "Flanking.Bivalent.TSS.Enh")]) 
          }
          priors[i,s] <<- p
          priors[s,i] <<- p
        }
      }
    });
  }

  ## SET SNP PRIOR
  # iterate over each of the snp genes and set prior according to
  # gtex pre-calculate priors
  snp.genes <- ranges$snp.genes$SYMBOL
  for(g in snp.genes) {
    # already filtered for sentinel id, just check gene
    idx <- which(grepl(paste0(paste0(",",g,"$"),"|", 
                              paste0(",",g,","),"|",
                              paste0("^",g,","),"|",
                              paste0("^",g,"$")), gtex.eqtl$gene_id))[1]
    # did we find the gene? then set prior
    if(!is.na(idx)){
      priors[g,id] <- priors[id,g] <- (1-gtex.eqtl[idx,"lfdr"])
    }
  }
  
  ## SET GENE-GENE PRIOR
  # here we add the priors based on the string interaction network
  genes <- colnames(priors)[!grepl("^rs|^cg", colnames(priors))];

  # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), genes),STRING.DB)
  edges <- edgeL(STRING.SUB)
  nodes <- nodes(STRING.SUB)
  edge.priors <- gtex.gg.cors
  rnames <- rownames(gtex.gg.cors)

  temp <- lapply(names(edges), function(n) {
    l <- edges[[n]]$edges;
    if (length(l) > 0) {
      temp2 <- sapply(l, function(i) {
        m <- nodes[i]
        # set priors for those validated connections found in STRING
        # check whether we have a specific prior for this gene-gene connection
        # match genes either at the beginning with a trailing "_" or at the end with a 
        # beginning "_" (rownames are: gene1_gene2)
        idxs <- which((grepl(paste0("^", n, "_"), rnames) & grepl(paste0("_", m, "$"), rnames)) |
                        (grepl(paste0("^", m, "_"), rnames) & grepl(paste0("_", n, "$"), rnames)))
        p <- gtex.gg.cors[idxs[1], "prior"]
        return(list(g1=m, g2=n, prior=p));
      });
      temp2
    }
  });

  temp <- unlist(temp)
  # did we have at least one connection?
  if(!is.null(temp)) {
    temp <- matrix(temp,ncol=3, byrow = T)
    for(i in 1:nrow(temp)) {
      g1 <- temp[i,1]
      g2 <- temp[i,2]
      prior <- as.numeric(temp[i,3])
      priors[g1,g2] <- prior
      priors[g2,g1] <- prior
    }
  }
  
  ## SET TF-CPG PRIOR
  # we will also set tf2cpg priors at this point, so get all TFs
  tfs <- ranges$tfs$SYMBOL
  # also get the chipseq context for our cpgs
  context <- get.chipseq.context(names(ranges$cpgs))
  
  # for all cpgs
  for(c in rownames(context)){
    for(tf in tfs) {
      # might be that the TF was measured in more than one cell line
      if(any(context[c,grepl(tf, colnames(context))])) {
        if(c %in% colnames(priors) &
           tf %in% colnames(priors)){
            priors[c,tf] <- 0.7
            priors[tf,c] <- 0.7
        } else {
          warning("TF-CpG link not available:",c,"-", tf, "\n")
        }
      }
    }
  }
  
  # define weights for our priors
  # do we need those?
  tf.cpg.weight <- 0.34
  cg.gene.weight <- 0.15;
  snp.gene.weight <- 0.15;
  gene.gene.weight <- 0.34;
  rest.weight <- (1 - (cg.gene.weight + snp.gene.weight + 
                       gene.gene.weight + tf.cpg.weight))

  cat("Prior weights used: \n")
  cat("\ttf.cpg:", tf.cpg.weight, "\n")
  cat("\tcg.gene:", cg.gene.weight, "\n")
  cat("\tsnp.gene:", snp.gene.weight, "\n")
  cat("\tgene.gene:", gene.gene.weight, "\n")
  cat("\trest.gene:", rest.weight, "\n")

  # build weighted prior matrix
  # TODO we surely could do this in a more refined way...
  priors.w <- priors
  for(j in 2:ncol(priors)){
    c <- colnames(priors)[j]
    for(i in 1:(nrow(priors)-1)){
      r <- rownames(priors)[i]
      
      # check which prior to use for this current link
      
      # cg-tf
      if(grepl("cg", c) & r %in% ranges$tfs$SYMBOL) {
        priors.w[i,j] <- priors[i,j] * tf.cpg.weight
      } else if(grepl("cg", c) & r %in% ranges$cpg.genes$SYMBOL) {
      # cg-gene
        priors.w[i,j] <- priors[i,j] * cg.gene.weight
      } else if(grepl("rs", c) & r %in% ranges$snp.genes$SYMBOL) { 
      # snp-gene
        priors.w[i,j] <- priors[i,j] * snp.gene.weight
      } else if(c %in% genes & r %in% genes) {
      # gene-gene
        priors.w[i,j] <- priors[i,j] * gene.gene.weight
      } else {
      # rest
        priors.w[i,j] <- priors[i,j] * rest.weight
      }
    }
  }
  colnames(priors.w) <- rownames(priors.w) <- colnames(priors)
  
  # copy upper tri matrix to lower tri
  for(i in 1:nrow(priors.w)) {
    for(j in 1:(i-1)) {
      priors.w[i,j] <- priors.w[j,i];
    }
  }

  # scale to overall increase priors.
  # TODO check whether this is even allowed. this would
  # increase the relevance of priors for the bdgraph algorithm?
  # Also: we found the influence of priors on model almost too large in 
  # a first run, consider this when thinking about uncommenting line below
  #priors.w <- priors.w / (max(priors.w) + 0.1)

  return(priors)
}

#' Gets priors for the GGMs from the GTEX dataset
#'
#' Creates binned priors based on the distances of SNPs and their cis genes as well
#' as for gene-gene connections. I.e. for the distance priors a window is binned
#' and for each bin the percentage of significant cis-associations to genes is
#' recorded. Since we want to transfer those to our data, we use the
#' currently set window size of our main script for this analysis.
#'
#' @param nbins The number of bins to create
#' @param p.cutoff The pvalue cutoff for significant associations
#' @param gtex.eqtl.file The data file in which to find the processed gtex associations
#' @param gtex.rpkm.file The file containing all rpkm values per GTEX sample and all genes
#' @param gtex.sampleDD.file The file containing the sample description for the samples
#' in the gtex.rpkm.file
#' @param chunk.size The number of lines to be read in one chunk (we expect quite a
#' huge file...)
#'
#' @return nothing
#'
#'
create.priors <- function(nbins, window) {

  source("R/lib.R")
  library(qvalue)
  library(data.table)
  library(graph)
  library(fdrtool)
  
  load.string.db();
  STRING.NODES <- nodes(STRING.DB)
  STRING.EDGES <- graph::edges(STRING.DB)

  gtex.eqtl.file <- paste0("data/current/gtex/Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt.gz")
  gtex.rpkm.file <- paste0("data/current/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz")
  gtex.sampleDS.file <- paste0("data/current/gtex/GTEx_Data_V6_Annotations_SampleAttributesDS.txt")
  gtex.phenotypeDS.file <- paste0("data/current/gtex/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt")
  
  if(!file.exists(gtex.eqtl.file)){
    stop("Provided gtex-eqtl file does not exist.\n");
  }
  if(!file.exists(gtex.rpkm.file)){
    stop("Provided gtex-rpkm file does not exist.\n");
  }
  if(!file.exists(gtex.sampleDS.file)){
    stop("Provided gtex-sample file does not exist.\n");
  }
  if(!file.exists(gtex.phenotypeDS.file)){
    stop("Provided gtex-sample file does not exist.\n");
  }

  GTEX.PLOTS <- paste0("results/current/gtex_plots/")
  dir.create(GTEX.PLOTS)
  
  # create the binned vector
  step <- window/nbins;
  bins <- seq(0,window-step, by=step);

  # read the gzfile chunk-wise
  # open handle
  f = paste0("zcat ", gtex.eqtl.file)

  # format: gene_id variant_id      tss_distance    pval_nominal    slope   slope_se
  pairs <- fread(f,
                 header=T, 
                 select=c("variant_id", "tss_distance", "pval_nominal"), 
                 data.table=F)
#  pairs <- cbind.data.frame(pairs, lfdr=qvalue(pairs$pval_nominal)$lfdr)
  
  # read the annotation information for the snps
#  annot <- fread("data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt",
#                 header=T,
#                 select=c("VariantID", "RS_ID_dbSNP135_original_VCF", "RS_ID_dbSNP142_CHG37p13"),
#                 data.table=F)
  
  # check whether we have annotation for all SNPs in the eQTL data
  # this seems not true, however there might be SNPs without rsID?
  # also: Check the extra column of rsIDs for overlap. maybe we can get all the
  # information if we use both columns?
  
 # all(pairs$variant_id %in% annot$VariantID)
  
  
  
  eqtl2bin <- findInterval(abs(pairs$tss_distance), bins)
  
  # for each bin, get the pvalues
  temp <- lapply(1:length(bins), function(b) {
    
    pvals <- as.numeric(pairs[which(eqtl2bin == b), "pval_nominal"])
    
    pdf(paste0(GTEX.PLOTS, "/gtex_bin", b, ".pval.pdf"))
    hist(pvals, breaks=100, main=paste0("bin ", b))
    dev.off()
    n <- 1-pi0est(pvals)$pi0
    return(n)
  })

  # relative amount of significant associations as priors
  gtex.eqtl.priors <- unlist(temp);
  # scale between 1 an 0
  gtex.eqtl.priors <- gtex.eqtl.priors / max(gtex.eqtl.priors);
  
  save(file=paste0("results/current/",
                   "/gtex.eqtl.priors.RData"), gtex.eqtl.priors)
  
  # plot the prior distribution (significant/bin)
  pdf(paste0(GTEX.PLOTS, "/gtex.eqtl.priors.pdf"))
  barplot(gtex.eqtl.priors, xlab="bins")
  dev.off()
  
  cat("Finished eQTL priors. Starting gene-gene priors now.\n");

  # format: 1:ID, 7:SMTSD (tissue)
  samples <- fread(gtex.sampleDS.file, select=c("SAMPID", 
                                                "SMTSD", 
                                                "SMNABTCH",
                                                "SMGEBTCH", 
                                                "SMRIN"));
  samples <- samples[which(samples$SMTSD=="Whole Blood"),]
  samples <- cbind.data.frame(id=samples$SAMPID, 
                              type=samples$SMTSD, 
                              batch1=samples$SMNABTCH, 
                              batch2=samples$SMGEBTCH,
                              RIN=samples$SMRIN,stringsAsFactors=F);
  # in order to be able to merge samples with covariate frame
  samples$donor_id <- gsub("(.{4}-.{4,5})-.*", "\\1", samples$id)
  # drop low-RIN samples and where RIN is NA 
  # see also: 
  # http://www.gtexportal.org/home/documentationPage#staticTextSampleQuality
  samples <- samples[complete.cases(samples),]
  samples <- samples[which(samples$RIN>=6),]
  # get covariates (age, sex) for all samples used
  covariates <- read.table(gtex.phenotypeDS.file,
                           header=T, 
                           sep="\t")
  colnames(covariates) <- c("id", "sex", "age", "death_hardy")
  
  samples <- merge(samples, covariates, by.x="donor_id", by.y="id")
  rownames(samples) <- samples$id
  
  cat("Processing", nrow(samples), "GTEX samples.\n")
  
  # now calculate the data-driven priors for gene-gene interactions
  rpkms <- fread(paste0("zcat ", gtex.rpkm.file), skip=2, sep="\t", 
                 header=T, select=c("Name","Description",samples$id));
  gene.matrix <- as.matrix(rpkms[,-c(1,2),with=F])
  rownames(gene.matrix) <- rpkms$Description
  rm(rpkms)
  
  cat("Subsetting genes.\n");

  gene.matrix <- t(gene.matrix[which(rownames(gene.matrix) %in% STRING.NODES),]);
  gene.matrix <- log10(gene.matrix+1)
  gene.matrix <- get.residuals(as.data.frame(gene.matrix),
                               samples[rownames(gene.matrix),])
  # plot gene expression values
  pdf(paste0(GTEX.PLOTS, "/expression.values.pdf"))
  hist(gene.matrix, breaks=150, xlab="expression residuals")
  dev.off()
  
  # calculate correlations for priors
  cat("Calculating gene-wise correlations.\n");
  
  cnames <- colnames(gene.matrix)

  temp <- lapply(STRING.NODES, function(e) {
    if(e %in% names(STRING.EDGES)) {
      temp2 <- lapply(STRING.EDGES[[e]], function(j) {
        if((e %in% cnames) & (j %in% cnames)) {
          x <- gene.matrix[,e];
          y <- gene.matrix[,j];
          # check sd
          if(sd(x)!=0 | sd(y) != 0) {
            pval <- cor.test(x,y)$p.value;
            rho <- cor(x,y)
            return(c(e,j,pval, rho))
          }
        }
      });
      return(temp2)
    }
  })

  cat("Done.\n");

  gtex.gg.cors <- matrix(unlist(temp), ncol=4, byrow = T)
  colnames(gtex.gg.cors) <- c("g1", "g2", "pval", "correlation")
  gtex.gg.cors <- as.data.frame(gtex.gg.cors, stringsAsFactors=F)  
  gtex.gg.cors$pval <- as.numeric(gtex.gg.cors$pval)
  gtex.gg.cors$correlation <- as.numeric(gtex.gg.cors$correlation)
  gtex.gg.cors$lfdr <- fdrtool(gtex.gg.cors$pval, statistic="pvalue")$lfdr
  gtex.gg.cors$prior <- 1-gtex.gg.cors$lfdr
  
  # report the pi0 and the proportion of significant tests for qvalue packge
  pi1 <- 1-pi0est(gtex.gg.cors[,"pval"])$pi0
  print(paste0("1-pi0 is: ", pi1))

  # report prior
  save(file=paste0("results/current/gtex.gg.cors.RData"), gtex.gg.cors)
}

#' Load gtex priors into environment
#'
#' @param sentinel The sentinel id. If given, eqtl data will immediately
#' be filtered for this SNP id. default: NULL
#' 
#' @author Johann Hawe
#' @return nothing
#'
load.gtex.priors <- function(sentinel=NULL) {
  gtex.eqtl <- readRDS("results/current/gtex.eqtl.priors.rds");
  load("results/current/gtex.gg.cors.RData");
  
  # keep only some of the columns of the eqtl priors
  gtex.eqtl <- gtex.eqtl[,c("gene_id",
                              "RS_ID_dbSNP135_original_VCF",
                              "RS_ID_dbSNP142_CHG37p13",
                              "lfdr")]
  
  # check whether to filter eqtl data
  
  if(!is.null(sentinel)){
    gtex.eqtl <- gtex.eqtl[(gtex.eqtl$RS_ID_dbSNP135_original_VCF==sentinel | 
                              gtex.eqtl$RS_ID_dbSNP142_CHG37p13==sentinel), ]
  }
  assign("gtex.eqtl", gtex.eqtl, .GlobalEnv);
  assign("gtex.gg.cors", gtex.gg.cors, .GlobalEnv);
}

#' Preprocesses the gtex Whole_Blood eqtl data to add a column of local FDR
#' values which we will use (1-lfdr) for our snp-gene priors
#' Output is written to a predetermined gz-file
#'
#' @author Johann Hawe
#' @return nothing
#'
preprocess.gtex.eqtl <- function(){
  library(fdrtool)
  library(data.table)
  library(Homo.sapiens)

  # load all gtex eqtl results
  gtex.eqtl.file <- paste0("data/current/gtex/Whole_Blood_Analysis.v6p.all_snpgene_pairs.txt.gz")
  # open handle
  f = paste0("zcat ", gtex.eqtl.file)
  
  # format: gene_id variant_id      tss_distance    pval_nominal    slope   slope_se
  pairs <- fread(f,
                 header=T, 
                 select=c("gene_id", "variant_id", "pval_nominal"), 
                 data.table=T)
  # get local fdr for all entries
  pairs[, ("lfdr") := fdrtool(pairs$pval_nominal, statistic="pvalue")$lfdr]

  # load annotation
  annot <- fread("data/current/gtex/GTEx_Analysis_v6_OMNI_genot_1KG_imputed_var_chr1to22_info4_maf01_CR95_CHR_POSb37_ID_REF_ALT.txt",
                 header=T,
                 select=c("VariantID", "RS_ID_dbSNP135_original_VCF", "RS_ID_dbSNP142_CHG37p13"),
                 data.table=T,
                 key="VariantID")
  
  # set the rsIDs for the pairs
  pairs[, ("RS_ID_dbSNP135_original_VCF") := annot[pairs$variant_id,mult="first"]$RS_ID_dbSNP135_original_VCF]
  pairs[, ("RS_ID_dbSNP142_CHG37p13") := annot[pairs$variant_id, mult="first"]$RS_ID_dbSNP142_CHG37p13]

  # get gene symbols
  ids <- sapply(strsplit(pairs$gene_id, "." ,fixed=T) ,"[", 1)
  symbols <- select(Homo.sapiens,
                     keys=unique(ids),
                     keytype="ENSEMBL",
                     columns="SYMBOL")
  symbols <- symbols[!is.na(symbols$SYMBOL),]
  
  pairs <- pairs[which(ids %in% symbols$ENSEMBL)]
  ids <- ids[ids %in% symbols$ENSEMBL]
  
  symbols.by.id <- tapply(symbols$SYMBOL, symbols$ENSEMBL, function(x) paste(x,collapse=","))
  pairs[, ("gene_id") := symbols.by.id[ids]]

  pairs[,"gene_id"] <- paste(symbols[ids,"SYMBOL"], collapse=",")

  # output gzip file
  out <- "results/current/gtex.eqtl.priors.rds"
  saveRDS(file=out, pairs)
}
