# Get priors for links between entities
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
  message("WARNING: Assuming single SNP in given data.")
  id <- nodes[grepl("^rs", nodes)]
  
  # load string db
  load.string.db();
  # load predefined prior definitions from gtex
  load.gtex.priors();
  
  # subset eqtl priors
  if(id %in% gtex.eqtl$RS_ID_dbSNP135_original_VCF | 
     id %in% gtex.eqtl$RS_ID_dbSNP142_CHG37p13) {
    gtex.eqtl <- gtex.eqtl[(gtex.eqtl$RS_ID_dbSNP135_original_VCF==id | 
                              gtex.eqtl$RS_ID_dbSNP142_CHG37p13==id), ]
  } else {
    cat("WARNING: Sentinel", id, "has no GTEx eQTL!\n")
    gtex.eqtl <- NULL
  }
  
  # default window (needed for distance calculation) and number of bins
  wnd <- 1e6;
  nbins <- 200;
   # create bin vector
  step <- wnd/nbins;
  bins <- seq(0,wnd-step, by=step);
  
  
  # get distances to all genes for cpgs  and snps
  cpg.dists <- get.nearby.ranges(ranges$cpgs, 
                                 promoters(ranges$cpg.genes))

  # sanity check
  if(!all(names(cpg.dists) %in% nodes)){
    cat("Some elements were not available in data:\n")
    na <- names(cpg.dists)
    cat(na[which(!na %in% nodes)])
    cpg.dists <- cpg.dists[!names(cpg.dists) %in% na]
    warning("WARNING: Some data mysteriously went missing..")
  }
  
  # now build the prior matrix using a pseudo prior
  pseudo.prior <- 1e-7
  priors <- matrix(data = pseudo.prior, nrow = length(nodes), ncol=length(nodes))
  colnames(priors) <- rownames(priors) <- nodes

  # load annotation needed for cpg2gene priors
  epigen.states <- read.table("data/current/epigenetic_state_annotation_weighted_all_sentinels.txt", 
                        header = T, sep="\t", row.names = 1)
  
  ## SET CPG-CPGGENE PRIOR
  for(i in names(cpg.dists)){
    gene.ranges <- cpg.dists[[i]];
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
    })
  }

  ## SET SNP PRIOR (if available)
  # iterate over each of the snp genes and set prior according to
  # gtex pre-calculated priors
  if(exists("gtex.eqtl")) {
    snp.genes <- ranges$snp.genes$SYMBOL
    for(g in snp.genes) {
      # already filtered for sentinel id, just check gene
      idx <- which(grepl(paste0(paste0(",",g,"$"),"|", 
                                paste0(",",g,","),"|",
                                paste0("^",g,","),"|",
                                paste0("^",g,"$")), gtex.eqtl$symbol))[1]
      # did we find the gene? then set prior
      if(!is.na(idx)){
        p <- (1-gtex.eqtl[idx]$lfdr)
        if(p <= pseudo.prior) {
          p <- pseudo.prior
        }
        if(p==1) { 
          p = 1 - pseudo.prior
        }
        priors[g,id] <- priors[id,g] <- p
      }
    }
  }
  
  ## SET GENE-GENE PRIOR
  # here we add the priors based on the string interaction network
  genes <- colnames(priors)[!grepl("^rs|^cg", colnames(priors))];

  # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), genes),STRING.DB)
  sedges <- edgeL(STRING.SUB)
  snodes <- nodes(STRING.SUB)
  edge.priors <- gtex.gg.cors
  rnames <- rownames(gtex.gg.cors)

  temp <- lapply(names(sedges), function(n) {
    l <- sedges[[n]]$edges;
    if (length(l) > 0) {
      temp2 <- sapply(l, function(i) {
        m <- snodes[i]
        # set priors for those validated connections found in STRING
        # check whether we have a specific prior for this gene-gene connection
        # match genes either at the beginning with a trailing "_" or at the end with a 
        # beginning "_" (rownames are: gene1_gene2)
        idxs <- which((grepl(paste0("^", n, "_"), rnames) & grepl(paste0("_", m, "$"), rnames)) |
                        (grepl(paste0("^", m, "_"), rnames) & grepl(paste0("_", n, "$"), rnames)))[1]
        if(!is.na(idxs)) {
          p <- pseudo.prior + gtex.gg.cors[idxs, "prior"]
        } else {
          p <- pseudo.prior
        }
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
      priors[g1,g2] <- priors[g2,g1] <- pseudo.prior + prior
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
          message("WARNING: TF-CpG link not available:",c,"-", tf, "\n")
        }
      }
    }
  }
  
  # define weights for our priors
  # do we need those?
  if(FALSE) {
    message("Using priors weighting.")
    tf.cpg.weight <- 0.3
    cg.gene.weight <- 0.2;
    snp.gene.weight <- 0.3;
    gene.gene.weight <- 0.19;
    rest.weight <- (1 - (cg.gene.weight + snp.gene.weight + 
                         gene.gene.weight + tf.cpg.weight))
  
    cat("Prior weights used: \n")
    cat("\ttf.cpg:", tf.cpg.weight, "\n")
    cat("\tcg.gene:", cg.gene.weight, "\n")
    cat("\tsnp.gene:", snp.gene.weight, "\n")
    cat("\tgene.gene:", gene.gene.weight, "\n")
    cat("\trest:", rest.weight, "\n")
  
    # build weighted prior matrix
    # TODO we surely could do this in a more refined way...
    priors.w <- priors
    for(j in 2:ncol(priors)){
      c <- colnames(priors)[j]
      for(i in 1:(nrow(priors)-1)){
        r <- rownames(priors)[i]
        
        p <- priors[i,j]
        
        # perform weighting only on set priors, skip  
        # pseudo prior values
        if(p != pseudo.prior) {
          # check which prior to use for this current link
          
          # cg-tf
          if(grepl("cg", c) & r %in% ranges$tfs$SYMBOL) {
            p <- p * tf.cpg.weight
          } else if(grepl("cg", c) & r %in% ranges$cpg.genes$SYMBOL) {
          # cg-gene
            p <- p * cg.gene.weight
          } else if(grepl("rs", c) & r %in% ranges$snp.genes$SYMBOL) { 
          # snp-gene
            p <- p * snp.gene.weight
          } else if(c %in% genes & r %in% genes) {
          # gene-gene
            p <- p * gene.gene.weight
          }
          priors.w[j,i] <- priors.w[i,j] <- max(pseudo.prior, p)
        }
      }
    }
    colnames(priors.w) <- rownames(priors.w) <- colnames(priors)
  
    # scale to overall increase priors.
    # TODO check whether this is even allowed. this would
    # increase the relevance of priors for the bdgraph algorithm?
    # Also: we found the influence of priors on model almost too large in 
    # a first run, consider this when thinking about uncommenting line below
    #priors.w <- priors.w / (max(priors.w) + pseudo.prior)
  } else {
    priors.w <- priors
  }
  # sanity check, matrix should not contain 0s or 1s or values smaller than our
  # peudo prior
  if(any(priors.w==0) | any(priors.w==1) | any(priors.w<pseudo.prior)) {
    stop(paste0("ERROR: Sanity check for priors failed. ",
                "(contain 0s/1s or values smaller than ", pseudo.prior, ")"))
  }

  return(priors.w)
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
                 header=T, select=c("Name","Description",samples$id))
  gene.matrix <- as.matrix(rpkms[,-c(1,2),with=F])
  rownames(gene.matrix) <- rpkms$Description
  rm(rpkms)
  
  cat("Subsetting genes.\n")

  gene.matrix <- t(gene.matrix)
  gene.matrix <- log2(gene.matrix+1)
  gene.matrix <- normalize.expression(gene.matrix)
  gene.matrix <- gene.matrix[,which(colnames(gene.matrix) %in% STRING.NODES)]
  # plot gene expression values
  pdf(paste0(GTEX.PLOTS, "/expression.values.pdf"))
  hist(gene.matrix, breaks=150, xlab="expression residuals")
  dev.off()
  
  # calculate correlations for priors
  cat("Calculating gene-wise correlations.\n")
  
  cnames <- colnames(gene.matrix)

  temp <- lapply(STRING.NODES, function(e) {
    if(e %in% names(STRING.EDGES)) {
      temp2 <- lapply(STRING.EDGES[[e]], function(j) {
        if((e %in% cnames) & (j %in% cnames)) {
          x <- gene.matrix[,e]
          y <- gene.matrix[,j]
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

  cat("Done.\n")

  gtex.gg.cors <- matrix(unlist(temp), ncol=4, byrow = T)
  colnames(gtex.gg.cors) <- c("g1", "g2", "pval", "correlation")
  rownames(gtex.gg.cors) <- with(gtex.gg.cors, paste(g1,g2,sep="_"))
  gtex.gg.cors <- as.data.frame(gtex.gg.cors, stringsAsFactors=F)  
  gtex.gg.cors$pval <- as.numeric(gtex.gg.cors$pval)
  gtex.gg.cors$correlation <- as.numeric(gtex.gg.cors$correlation)
  gtex.gg.cors$lfdr <- fdrtool(gtex.gg.cors$pval, statistic="pvalue")$lfdr
  gtex.gg.cors$prior <- 1-gtex.gg.cors$lfdr
  
  # plot correlation histogram
  pdf(paste0(GTEX.PLOTS, "/expression.correlation.pdf"))
  hist(gtex.gg.cors$correlation, breaks=100, xlab="correlation", main="gtex gene corerlations")
  abline(v=0, col="red")
  dev.off()
  
  # report the pi0 and the proportion of significant tests for qvalue packge
  pi1 <- 1-pi0est(gtex.gg.cors[,"pval"])$pi0
  print(paste0("1-pi0 is: ", pi1))

  # report prior
  saveRDS(file=paste0("results/current/gtex.gg.cors.rds"), gtex.gg.cors)
}

#' Load gtex priors into environment
#'
#' @param sentinel The sentinel id. If given, eqtl data will immediately
#' be filtered for this SNP id. default: NULL
#' 
#' @author Johann Hawe
#' @return nothing
#'
load.gtex.priors <- function(sentinel=NULL, env=.GlobalEnv) {
  # check whether prior were already loaded
  if(with(env, exists("gtex.eqtl")) &
     with(env, exists("gtex.gg.cors"))){
    # nothing to do
    return();
  }
  cat("Loading gtex priors.\n")
  gtex.eqtl <- readRDS("results/current/gtex.eqtl.priors.rds");
  gtex.gg.cors <- readRDS("results/current/gtex.gg.cors.rds");
  
  # keep only some of the columns of the eqtl priors
  gtex.eqtl <- gtex.eqtl[,c("symbol",
                              "RS_ID_dbSNP135_original_VCF",
                              "RS_ID_dbSNP142_CHG37p13",
                              "lfdr")]
  
  # check whether to filter eqtl data
  
  if(!is.null(sentinel)){
    if(sentinel %in% gtex.eqtl$RS_ID_dbSNP135_original_VCF | 
       sentinel %in% gtex.eqtl$RS_ID_dbSNP142_CHG37p13) {
      gtex.eqtl <- gtex.eqtl[(gtex.eqtl$RS_ID_dbSNP135_original_VCF==sentinel | 
                              gtex.eqtl$RS_ID_dbSNP142_CHG37p13==sentinel), ]
    } else {
      cat("WARNING: Sentinel", sentinel, "has no GTEx eQTL!\n")
      gtex.eqtl <- NULL
    }
  }
  # check whether we have eqtl results
  if(!is.null(gtex.eqtl)){
    assign("gtex.eqtl", gtex.eqtl, env);
  }
  assign("gtex.gg.cors", gtex.gg.cors, env);
}

#' Preprocesses the gtex Whole_Blood eqtl data to add a column of local FDR
#' values which we will use (1-lfdr) for our snp-gene priors
#' Output is written to a predetermined gz-file
#'
#' @author Johann Hawe
#' @return nothing
#'
preprocess.gtex.eqtl <- function() {
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
  pairs[, ("symbol") := symbols.by.id[ids]]

  # output gzip file
  out <- "results/current/gtex.eqtl.priors.rds"
  saveRDS(file=out, pairs)
}

#' Method creates priors for TF~CpG links for application of the GGM
#' Loads Remap/ENCODE TFBS and uses score column of bed files to infer 
#' prior values.
#' 
#' TODO: a possible extension could be to get the sequences where the TFBS are 
#' identified to be bound and check their PWM score for those positions...
#' 
#' @author Johann Hawe
#' 
create.tfbs.priors <- function() {
  ## annotation of CpG sites with functional genomics data and TFBS predictions
  
  ## compute TF affinities with the tRap package
  library(tRap)
  
  ## load the CpG positions from the bioconductor package
  if (!require(FDb.InfiniumMethylation.hg19)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("FDb.InfiniumMethylation.hg19")
  }
  
  ## genome sequences
  if (!require(BSgenome.Hsapiens.UCSC.hg19)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
  }
  
  ## get 100bp intervals around the CpG sites
  cpgs = features(FDb.InfiniumMethylation.hg19)
  context = resize(cpgs, 100, fix="center")
  
  ## get the corresponding sequence
  seq = as.character(getSeq(Hsapiens, context))
  
  ## Code below would presumably performs the PWM steps, however
  ## currently only for a single TF
  ## get the PWMs
  #data(transfac)
  #matrices = read.transfac("data/current/pwms/znf333.txt")
  ## onlly use vertebrate factors
  #matrices = matrices[substr(names(matrices), 1, 1) == "V"]
  #affinities = mclapply(matrices, function(pwm) affinity(pwm, seq, both.strands=T))
  #save(affinities, file="results/current/cpgs_with_affinties.RData")
  
  
  ## also annotate CpGs with ChIP-seq binding sites
  library(rtracklayer)
  library(data.table)
  
  tfbs = import("data/current/tfbs/filPeaks_public.bed")
  ann = t(matrix(unlist(strsplit(values(tfbs)[,"name"], ".", fixed=T)), nrow=3))
  ann <- cbind(ann, score=tfbs$score)
  colnames(ann) = c("geo_id", "TF", "condition", "score")
  values(tfbs) = DataFrame(name=values(tfbs)[,"name"], data.frame(ann, stringsAsFactors=F))
  
  ## we write out a table with all conditions and select the blood related ones
  conditions = t(matrix(unlist(strsplit(unique(values(tfbs)[,"name"]), ".", fixed=T)), nrow=3))
  colnames(conditions) = c("geo_id", "TF", "condition")
  conditions = conditions[order(conditions[,"condition"]),]
  conditions = conditions[,c(1,3)]
  conditions = conditions[!duplicated(paste(conditions[,1], conditions[,2])),]
  conditions = data.frame(conditions, blood.related=F)
  # set the 'blood.related' flag
  for (term in c("amlpz12_leukemic", "aplpz74_leukemia", "bcell", "bjab", "bl41", "blood", "lcl", "erythroid", "gm", "hbp", "k562", "kasumi", "lymphoblastoid", "mm1s", "p493", "plasma", "sem", "thp1", "u937")) {
    conditions[grep(term, conditions[,2]),"blood.related"] = TRUE
  }
  
  # gets us 35 distinct TFs only!
  selected = tfbs[values(tfbs)[,"condition"] %in% conditions[conditions[,"blood.related"],"condition"]]
  selected$score <- as.numeric(selected$score)
  
  # gets us 24 distinct TFs only!
  selected <- selected[selected$score>0]
 
  # TODO!!
  
   ## load the encode tfs separately
  #encode = as.data.frame(fread("data/current/tfbs/wgEncodeRegTfbsClusteredWithCellsV3.bed", header=F))
  #encode = GRanges(seqnames=encode[,1], 
  #                 ranges=IRanges(encode[,2] + 1, 
  #                                encode[,3]), 
  #                 name=paste("ENCODE", encode[,4], tolower(encode[,6]), sep="."), 
  #                 geo_id="ENCODE", 
  #                 TF=encode[,4], 
  #                 condition=tolower(encode[,6]))
  
  #encode.lcl = encode[grep("gm", values(encode)[,"condition"])]
  #values(encode.lcl)[,"condition"] = "lcl"
  #encode.k562 = encode[grep("k562", values(encode)[,"condition"])]
  #values(encode.k562)[,"condition"] = "k562"
  #selected = c(selected, encode.lcl, encode.k562)
  #chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
  #chip.exp = unique(chip)
  
  ## create an annotation matrix for the CpGs
  #tfbs.ann = sapply(chip.exp, function(x) overlapsAny(context, selected[chip == x]))
  #rownames(tfbs.ann) = names(cpgs)
  
  #saveRDS(tfbs.ann, file="results/current/cpgs_with_chipseq_context_100.rds")
  
}
