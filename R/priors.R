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
  
  # load string db
  load.string.db();
  # load predefined prior definition
  load.gtex.priors();
  
  # assume single SNP
  id <- nodes[grepl("^rs", nodes)]
  
  # default window (needed for distance calculation) and number of bins
  wnd <- 1e6;
  nbins <- 200;
   # create bin vector
  step <- wnd/nbins;
  bins <- seq(0,wnd-step, by=step);
  
  
  # get distances to all genes for cpgs and snps
  cpg.dists <- get.nearby.ranges(ranges$cpgs, 
                                 promoters(ranges$cpg.genes))
  temp <- ranges$snp.genes
  temp$distance <- distance(ranges$sentinel.range,
                                 promoters(ranges$snp.genes));
  snp.dists <- list()
  snp.dists[[id]] <- temp
  
  all.dists <- append(cpg.dists, snp.dists);

  # sanity check
  if(!all(names(all.dists) %in% nodes)){
    cat("Some elements were not available in data:\n")
    na <- names(all.dists)
    cat(na[which(!na %in% nodes)])
    all.dists <- all.dists[!names(all.dists) %in% na]
    warning("Some data misteriously went missing..")
  }
  
  # now build the prior matrix using a pseudo prior
  pseudo.prior <- 0.0000001
  priors <- matrix(data = pseudo.prior, nrow = length(nodes), ncol=length(nodes))
  colnames(priors) <- rownames(priors) <- nodes

  # load annotation needed for cpg2gene priors
  epigen.states <- read.table("data/current/epigenetic_state_annotation_weighted_all_sentinels.txt", 
                        header = T, sep="\t", row.names = 1)
  
  for(i in names(all.dists)){

    gene.ranges <- all.dists[[i]];

    # indicator whether the current distance is related to a cpg
    isCpGDist <- grepl("^cg", i);
    
    temp <- lapply(gene.ranges, function(r) {
      s <- r$SYMBOL;
      if(s %in% colnames(priors)) {
        # set the basic prior based on distance (the larger the distance, the lower the prior should be)
        # but scale it to be between 0 and 1
        dist <- r$distance
        if(!is.na(dist)){
          # distance based priors for snps (gtex eQTLs)
          if(!isCpGDist) {
             idx <- findInterval(dist, bins);
             p <- gtex.eqtl.priors[idx];
          }
          else {
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
          }
  
          priors[i,s] <<- p
          priors[s,i] <<- p
        }
      }
    });
  }

  # here we add the priors based on the string interaction network
  genes <- colnames(priors)[!grepl("^rs|^cg", colnames(priors))];

  # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(STRING.DB), genes),STRING.DB)
  edges <- edgeL(STRING.SUB)
  nodes <- nodes(STRING.SUB)
  edges.prior <- gtex.gg.prior
  
  temp <- lapply(names(edges), function(n) {
    l <- edges[[n]]$edges;
    if (length(l) > 0) {
      temp2 <- sapply(l, function(i) {
        m <- nodes[i]
        # set priors for those validated connections found in STRING
        return(list(g1=m, g2=n, prior=edges.prior));
      });
      temp2
    }
  });

  temp <- unlist(temp)
  # did we have at least one connection?
  if(!is.null(temp)) {
    temp <- matrix(temp,ncol=3, byrow = T);

    for(i in 1:nrow(temp)) {
      g1 <- temp[i,1];
      g2 <- temp[i,2];
      prior <- as.numeric(temp[i,3])
      priors[g1,g2] <- prior;
      priors[g2,g1] <- prior;
    }
  }
  
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
  # do we need those? discard for now
  #tf.cpg.weight <- 0.34
  #cg.gene.weight <- 0.15;
  #snp.gene.weight <- 0.15;
  #gene.gene.weight <- 0.34;
  #rest.weight <- 1- (cg.gene.weight + snp.gene.weight + 
  #                     gene.gene.weight + tf.cpg.weight)

  #cat("Prior weights used: \n")
  #cat("\ttf.cpg:", tf.cpg.weight, "\n")
  #cat("\tcg.gene:", cg.gene.weight, "\n")
  #cat("\tsnp.gene:", snp.gene.weight, "\n")
  #cat("\tgene.gene:", gene.gene.weight, "\n")
  #cat("\trest.gene:", rest.weight, "\n")

  # build weighted prior matrix
  #priors.w <- priors
  #for(j in 2:ncol(priors)){
  #  c <- colnames(priors)[j]
  #  for(i in 1:(nrow(priors)-1)){
  #    r <- rownames(priors)[i]
      
      # check which prior to use for this current link
      
      # cg-tf
  #    if(grepl("cg", c) & r %in% sentinel$enriched.tfs$SYMBOL) {
  #      priors.w[i,j] <- priors[i,j] * tf.cpg.weight
  #    } else if(grepl("cg", c) & r %in% sentinel$cpg.genes$SYMBOL) {
      # cg-gene
  #      priors.w[i,j] <- priors[i,j] * cg.gene.weight
  #    } else if(grepl("rs", c) & r %in% sentinel$snp.genes$SYMBOL) { 
      # snp-gene
  #      priors.w[i,j] <- priors[i,j] * snp.gene.weight
  #    } else if(c %in% genes & r %in% genes) {
      # gene-gene
  #      priors.w[i,j] <- priors[i,j] * gene.gene.weight
  #    } else {
      # rest
  #      priors.w[i,j] <- priors[i,j] * rest.weight
  #    }
  #  }
  #}
  #colnames(priors.w) <- rownames(priors.w) <- colnames(priors)
  
  # copy upper tri matrix to lower tri
  #for(i in 1:nrow(priors.w)) {
  #  for(j in 1:(i-1)) {
  #    priors.w[i,j] <- priors.w[j,i];
  #  }
  #}

  # scale to overall increase priors.
  # TODO check whether this is even allowed. this would
  # increase the relevance of priors for the bdgraph algorithm?
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

  library(qvalue)
  library(data.table)

  load.string.db();
 
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
                 select=c("tss_distance", "pval_nominal"), 
                 data.table=F)

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
  gene.matrix <- log10(gene.matrix+10)
  gene.matrix <- get.residuals(as.data.frame(gene.matrix),
                               samples[rownames(gene.matrix),])
  
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
            return(c(e,j,pval))
          }
        }
      });
      return(temp2)
    }
  })

  cat("Done.\n");

  gtex.gg.cors <- matrix(unlist(temp),ncol=3,byrow = T)
  colnames(gtex.gg.cors) <- c("g1", "g2", "cor.test.pval");
  gtex.gg.cors <- as.data.frame(gtex.gg.cors, stringsAsFactors=F)  
  gtex.gg.cors$cor.test.pval <- as.numeric(gtex.gg.cors$cor.test.pval)

  # get the pi0 and the proportion of significant tests for qvalue packge
  pi1 <- 1-pi0est(gtex.gg.cors[,3])$pi0
  
  # relative amount of significant associations
  gtex.gg.prior <- pi1
  
  # report prior
  save(file=paste0("results/current/gtex.gg.prior.RData"), gtex.gg.cors,
       gtex.gg.prior)
}

#' Load gtex priors into environment
#'
#' @return nothing
#'
load.gtex.priors <- function() {
  load("results/current/gtex.eqtl.priors.RData");
  load("results/current/gtex.gg.prior.RData");
  
  assign("gtex.eqtl.priors", gtex.eqtl.priors, .GlobalEnv);
  assign("gtex.gg.prior", gtex.gg.prior, .GlobalEnv);
}

