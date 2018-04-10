#' Get priors for links between entities ---------------------------------------
#'
#' Given appropriate information, creates a matrix for all GGM entitites (i.e.
#' snps, cpgs and cis-/trans-genes) reflecting respective prior probabilities 
#' for links between them (currently only based on 
#' distances of snps/cpgs <-> genes)
#'
#' @param ranges Set of ranges with the annotated data (e.g. snp.ranges, 
#' cpgs etc.) for which to extract the link priors
#'
#' @param nodes All the nodes/entities being analyzed.
#' 
#' @return A square matrix of link-priors
#'------------------------------------------------------------------------------
get.link.priors <- function(ranges, nodes, string_db) {
  
  # assume single SNP
  message("WARNING: Assuming single SNP in given data.")
  id <- nodes[grepl("^rs", nodes)]
  
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
  
  # now build the prior matrix using a pseudo prior
  pseudo_prior <- 1e-7
  priors <- matrix(data = pseudo_prior, 
                   nrow = length(nodes), 
                   ncol=length(nodes))
  colnames(priors) <- rownames(priors) <- nodes
  
  # load annotation needed for cpg2gene priors
  # TODO add those states for cpg priors once we got our storage back
  #epigen.states <- read.table("data/current/epigenetic_state_annotation_weighted_all_sentinels.txt", 
  #                      header = T, sep="\t", row.names = 1)
  
  # TODO for new structure of ranges!
  ## SET CPG-CPGGENE PRIOR
  # for(i in names(cpg.dists)){
  #   gene.ranges <- cpg.dists[[i]];
  #   temp <- lapply(gene.ranges, function(r) {
  #     s <- r$SYMBOL;
  #     if(s %in% colnames(priors)) {
  #       # set the basic prior based on distance (the larger the distance, the lower the prior should be)
  #       # but scale it to be between 0 and 1
  #       dist <- r$distance
  #       if(!is.na(dist)){
  #         # if the cpg is in range of the tss (within 200bp), 
  #         # set a specific prior for active TSS... 
  #         p <- pseudo.prior
  #         if(dist <= 200){
  #           # set the cpg.state prior
  #           p <- p + sum(epigen.states[i, c("Active.TSS", 
  #                                       "Flanking.Active.TSS", 
  #                                       "Bivalent.Poised.TSS",
  #                                       "Flanking.Bivalent.TSS.Enh")]) 
  #         }
  #         priors[i,s] <<- p
  #         priors[s,i] <<- p
  #       }
  #     }
  #   })
  # }
  
  ## SET SNP PRIOR (if available)
  # iterate over each of the snp genes and set prior according to
  # gtex pre-calculated priors
  if(exists("gtex.eqtl")) {
    snp_genes <- ranges$snp_genes$SYMBOL
    for(g in snp_genes) {
      # already filtered for sentinel id, just check gene
      idx <- which(grepl(paste0(paste0(",",g,"$"),"|", 
                                paste0(",",g,","),"|",
                                paste0("^",g,","),"|",
                                paste0("^",g,"$")), gtex.eqtl$symbol))[1]
      # did we find the gene? then set prior
      if(!is.na(idx)){
        p <- (1-gtex.eqtl[idx]$lfdr)
        if(p <= pseudo_prior) {
          p <- pseudo_prior
        }
        if(p==1) { 
          p = 1 - pseudo_prior
        }
        priors[g,id] <- priors[id,g] <- p
      }
    }
  }
  
  ## SET GENE-GENE PRIOR
  # here we add the priors based on the string interaction network
  genes <- colnames(priors)[!grepl("^rs|^cg", colnames(priors))]
  
  # get subset of edges which are in our current graph
  STRING.SUB <- subGraph(intersect(nodes(string_db), genes), string_db)
  sedges <- edgeL(STRING.SUB)
  snodes <- nodes(STRING.SUB)
  rnames <- rownames(gtex.gg.cors)
  
  temp <- lapply(names(sedges), function(n) {
    l <- sedges[[n]]$edges;
    if (length(l) > 0) {
      temp2 <- sapply(l, function(i) {
        m <- snodes[i]
        # set priors for those validated connections found in STRING
        # check whether we have a specific prior for this gene-gene connection
        # match genes either at the beginning with a trailing "_" or
        # at the end with a beginning "_" (rownames are: gene1_gene2)
        idxs <- which((grepl(paste0("^", n, "_"), rnames) & 
                         grepl(paste0("_", m, "$"), rnames)) |
                        (grepl(paste0("^", m, "_"), rnames) & 
                           grepl(paste0("_", n, "$"), rnames)))[1]
        if(!is.na(idxs)) {
          p <- pseudo_prior + gtex.gg.cors[idxs, "prior"]
        } else {
          p <- pseudo_prior
        }
        return(list(g1=m, g2=n, prior=p))
      })
      temp2
    }
  })
  
  temp <- unlist(temp)
  # did we have at least one connection?
  if(!is.null(temp)) {
    temp <- matrix(temp,ncol=3, byrow = T)
    for(i in 1:nrow(temp)) {
      g1 <- temp[i,1]
      g2 <- temp[i,2]
      prior <- as.numeric(temp[i,3])
      priors[g1,g2] <- priors[g2,g1] <- pseudo_prior + prior
    }
  }
  
  ## SET TF-CPG PRIOR
  # we will also set tf2cpg priors at this point, so get all TFs
  # tfs <- ranges$tfs$SYMBOL
  # # also get the chipseq context for our cpgs
  # context <- get.chipseq.context(names(ranges$cpgs))
  # 
  # # for all cpgs
  # for(c in rownames(context)){
  #   for(tf in tfs) {
  #     # might be that the TF was measured in more than one cell line
  #     if(any(context[c,grepl(tf, colnames(context))])) {
  #       if(c %in% colnames(priors) &
  #          tf %in% colnames(priors)){
  #           priors[c,tf] <- 0.7
  #           priors[tf,c] <- 0.7
  #       } else {
  #         message("WARNING: TF-CpG link not available:",c,"-", tf, "\n")
  #       }
  #     }
  #   }
  # }


  # sanity check, matrix should not contain 0s or 1s or values smaller than our
  # peudo prior
  if(any(priors==0) | any(priors==1) | any(priors<pseudo_prior)) {
    stop(paste0("ERROR: Sanity check for priors failed. ",
                "(contain 0s/1s or values smaller than ", pseudo_prior, ")"))
  }
  
  return(priors)
}

#' Load gtex priors into environment--------------------------------------------
#'
#' @param sentinel The sentinel id. If given, eqtl data will immediately
#' be filtered for this SNP id. default: NULL
#' 
#' @author Johann Hawe
#' @return nothing
#'------------------------------------------------------------------------------
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

#' Create priors for the GGMs from the GTEX dataset ----------------------------
#'
#' @author Johann Hawe
#' 
#' @return nothing
#'------------------------------------------------------------------------------
create.priors <- function(feqtl, fsnpInfo, frpkm, 
                          fsampleDS, fphenotypeDS, dplots, string_db) {
  
  # might be the plotting directory does not exist yet
  dir.create(dplots)
  
  # create the gene-gene priors
  create.gtex.gene.priors(fsampleDS, fphenotypeDS, frpkm, dplots, string_db)
  
  gc()
  
  # create the snp-gene priors
  create.gtex.eqtl.priors(feqtl, fsnpInfo)
}

#' Creates the gtex gene-gene priors under consideration of the ----------------
#' STRING database.
#' 
#' @param sampleDS.file The gtex sample information file
#' @param pheno.file The gtex phenotype file
#' @param rpkm.file File containing the gtex rpkm values
#' @param plot.dir Directory where to save summary plots
#' 
#' @author Johann Hawe
#'------------------------------------------------------------------------------
create.gtex.gene.priors <- function(sampleDS.file, pheno.file, 
                                    rpkm.file, plot.dir, string_db) {
  
  # load the string information
  STRING.NODES <- nodes(string_db)
  STRING.EDGES <- graph::edges(string_db)
  
  # create the Gene-Gene priors
  # format: 1:ID, 7:SMTSD (tissue)
  samples <- fread(sampleDS.file, select=c("SAMPID", 
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
  covariates <- read.table(pheno.file,
                           header=T, 
                           sep="\t")
  colnames(covariates) <- c("id", "sex", "age", "death_hardy")
  
  samples <- merge(samples, covariates, by.x="donor_id", by.y="id")
  rownames(samples) <- samples$id
  
  print(paste0("Processing ", nrow(samples), " GTEX samples."))
  
  # now calculate the data-driven priors for gene-gene interactions
  rpkms <- fread(paste0("zcat ", rpkm.file), skip=2, sep="\t", 
                 header=T, select=c("Name","Description",samples$id))
  gene.matrix <- as.matrix(rpkms[,-c(1,2),with=F])
  rownames(gene.matrix) <- rpkms$Description
  
  print("Subsetting genes.")
  
  gene.matrix <- t(gene.matrix)
  gene.matrix <- log2(gene.matrix+1)
  gene.matrix <- normalize.expression(gene.matrix)
  gene.matrix <- gene.matrix[,which(colnames(gene.matrix) %in% STRING.NODES)]
  # plot gene expression values
  pdf(paste0(plot.dir, "/expression.values.pdf"))
  hist(gene.matrix, breaks=150, xlab="expression residuals")
  dev.off()
  
  print("Calculating gene-wise correlations.\n")
  
  cnames <- colnames(gene.matrix)
  
  temp <- lapply(STRING.NODES, function(e) {
    if(e %in% names(STRING.EDGES)) {
      temp2 <- lapply(STRING.EDGES[[e]], function(j) {
        if((e %in% cnames) & (j %in% cnames)) {
          x <- gene.matrix[,e]
          y <- gene.matrix[,j]
          # check sd
          if(sd(x)!=0 | sd(y) != 0) {
            pval <- cor.test(x,y)$p.value
            rho <- cor(x,y)
            return(c(e,j,pval, rho))
          }
        }
      })
      return(temp2)
    }
  })
  
  print("Done.")
  
  gtex.gg.cors <- matrix(unlist(temp), ncol=4, byrow = T)
  colnames(gtex.gg.cors) <- c("g1", "g2", "pval", "correlation")
  gtex.gg.cors <- as.data.frame(gtex.gg.cors, stringsAsFactors=F)  
  rownames(gtex.gg.cors) <- with(gtex.gg.cors, paste(g1,g2,sep="_"))
  gtex.gg.cors$pval <- as.numeric(gtex.gg.cors$pval)
  gtex.gg.cors$correlation <- as.numeric(gtex.gg.cors$correlation)
  gtex.gg.cors$lfdr <- fdrtool(gtex.gg.cors$pval, statistic="pvalue")$lfdr
  gtex.gg.cors$prior <- 1-gtex.gg.cors$lfdr
  
  # plot correlation histogram
  pdf(paste0(plot.dir, "/expression.correlation.pdf"))
  hist(gtex.gg.cors$correlation, breaks=100, xlab="correlation", main="gtex gene corerlations")
  abline(v=0, col="red")
  dev.off()
  
  # report the pi0 and the proportion of significant tests for qvalue packge
  pi1 <- 1-pi0est(gtex.gg.cors[,"pval"])$pi0
  print(paste0("1-pi0 is: ", pi1))
  
  # report prior
  saveRDS(file=paste0("results/current/gtex.gg.cors.rds"), gtex.gg.cors)
}

#' Preprocesses the gtex Whole_Blood eqtl data to add a column of local FDR
#' values which we will use (1-lfdr) for our snp-gene priors
#' Output is written to a predetermined gz-file
#'
#' @author Johann Hawe
#' @return nothing
#'
create.gtex.eqtl.priors <- function(eqtl.file, snpInfo.file) {
  
  # open handle
  f = paste0("zcat ", eqtl.file)
  
  # format: gene_id variant_id      tss_distance    pval_nominal    slope   slope_se
  pairs <- fread(f,
                 header=T, 
                 select=c("gene_id", "variant_id", "pval_nominal"), 
                 data.table=T)
  # get local fdr for all entries
  pairs[, ("lfdr") := fdrtool(pairs$pval_nominal, statistic="pvalue")$lfdr]
  
  # load annotation
  annot <- fread(snpInfo.file,
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
  
  # output rds file
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

randomize_priors <- function(priors, percent) {
  library(reshape2)
  
  # get simple list of priors
  p <- priors
  p[upper.tri(p,T)] <- NA
  p <- melt(p, na.rm=T)
  p <- p[sample(1:nrow(p)),]
  
  # scale all non-pseudo priors to be within [0..1]
  p$scaled <- NA
  idxs <- which(p$value>min(p$value))
  p[idxs,]$scaled <- p[idxs,"value"]/sum(p[idxs,"value"])
  
  # identify the edges to randomize
  # we try length(idxs) times to get it as close as possible
  # to our desired range
  # define epsilon range
  if(percent > 0) {
    eps <- 0.01
    lo <- percent - eps
    up <- percent + eps
    best <- c()
    best_diff <- 1
    for(i in 1:length(idxs)*2) {
      rand <- cumsum_in_range(p, sample(idxs), lo, up)
      if(!is.null(rand)) {
        s <- sum(rand$scaled)
        d <- abs(s-percent)
        if(d < best_diff) {
          best <- rand
          best_diff <- d
        }
        if(d==0) {
          break
        }
      }
    }
    if(is.null(best)) {
      warning("No prior-switchings selected, returning default priors.")
      return(priors)
    }
    # now randomize the indices
    idxs_noprior <- which(is.na(p$scaled))
    idxs_noprior <- sample(idxs_noprior, nrow(best))
    
    # create list of switches to be performed
    switch <- cbind(prior=best$idx, no_prior=idxs_noprior)
    # switch in our matrix
    pseudo <- min(priors)
    for(i in 1:nrow(switch)) {
      p1 <- p[switch[i, "prior"], "Var1"]
      p2 <- p[switch[i, "prior"], "Var2"]
      np1 <- p[switch[i, "no_prior"], "Var1"]
      np2 <- p[switch[i, "no_prior"], "Var2"]
      
      # get prior value
      v <- priors[p1,p2]
      
      # set pseudo prior connection to new prior value
      priors[np1,np2] <- priors[np2,np1] <- v
      # set prior connection to pseudo value
      priors[p1,p2] <- priors[p2,p1] <- pseudo
    }
  }
  return(priors)
}

cumsum_in_range <- function(df, idxs, lo, up, col_name="scaled") {

  # get the  center value for exact matching
  percent <- (lo+up)/2
  
  if(percent == 0) {
    return(NULL)
  }
  rand <- c()
  for(i in idxs) {
    pi <- df[i,,drop=F]
    pi <- cbind(pi, idx=i)
    v <- pi[,col_name]
    s <- sum(c(rand[,col_name],v))
    
    # did not reach lower bound yet
    if(s<lo) {
      rand <- rbind(rand,pi)
      next
    }
    # above lower bound
    if(s>=lo) {
      # match exactly? then done
      if(s==percent) {
        rand <- rbind(rand,pi)
        break
      }
      # below upper bound
      if(s<=up) {
        rand <- rbind(rand,pi)
        # first above exact? then done
        if(s>percent) {
          break
        }
      } else {
        # exceeded with current value, stop
        break
      }
    }
  }
  return(rand)
}
