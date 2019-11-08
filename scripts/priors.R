#'------------------------------------------------------------------------------
#' Get priors for links between entities
#'
#' Given appropriate information, creates a matrix for all GGM entitites (i.e.
#' snps, cpgs and cis-/trans-genes) reflecting respective prior probabilities
#' for links between them (currently only based on
#' distances of snps/cpgs <-> genes). This methods assumes that the gtex prior
#' information has already been loaded
#'
#' @param ranges Set of ranges with the annotated data (e.g. snp.ranges,
#' cpgs etc.) for which to extract the link priors
#' @param nodes All the nodes/entities being analyzed.
#' @param ppi_db The loaded PPI db graph to be used (graphNEL)
#' @param fcpgcontext The CpG TF annotation
#' @param fcpg_annot The epigenetic state annotation file for the CpGs
#'
#' @return A square matrix of link-priors
#'------------------------------------------------------------------------------
get_link_priors <- function(ranges, nodes, ppi_db, fcpgcontext, fcpg_annot) {

  # ----------------------------------------------------------------------------
  # prepare some needed data
  # ----------------------------------------------------------------------------

  # assume single SNP
  print("Assuming single SNP in given data.")
  id <- nodes[grepl("^rs", nodes)]

  # ----------------------------------------------------------------------------
  # now build the prior matrix, set pseudo prior
  # ----------------------------------------------------------------------------
  pseudo_prior <- 1e-7
  priors <- matrix(data = pseudo_prior,
                   nrow = length(nodes),
                   ncol=length(nodes))
  colnames(priors) <- rownames(priors) <- nodes

  # ----------------------------------------------------------------------------
  # check whether we need CpG related priors
  # ----------------------------------------------------------------------------
  if(ranges$seed == "meqtl") {
    print("Setting CpG-gene priors")
    # load annotation needed for cpg2gene priors
    epigen.states <- read.table(fcpg_annot,
                                header = T, sep="\t", row.names = 1)

    # --------------------------------------------------------------------------
    ## SET CPG-CPGGENE PRIOR
    # --------------------------------------------------------------------------
    for(n in names(ranges$cpg_genes_by_cpg)){
      cpg <- ranges$cpgs[n]
      cgs <- ranges$cpg_genes_by_cpg[[n]]

      for(i in 1:length(cgs)) {
        g <- cgs[i]
        symbol <- g$SYMBOL
        d <- distance(g, cpg)
        if(symbol %in% colnames(priors)) {
          # set the basic prior based on distance
          # (the larger the distance, the lower the prior should be)
          # but scale it to be between 0 and 1
          if(!is.na(d)){
            # if the cpg is in range of the tss (within 200bp),
            # set a specific prior for active TSS...
            p <- pseudo_prior
            if(d <= 200){
              # set the cpg.state prior
              p <- p + sum(epigen.states[n, c("Active.TSS",
                                                "Flanking.Active.TSS",
                                                "Bivalent.Poised.TSS",
                                                "Flanking.Bivalent.TSS.Enh")])
            }
            if(p>=1) p <- 1 - pseudo_prior
            priors[n,symbol] <- priors[symbol,n] <- p
          } else {
            stop("CpG-gene not proximal to CpG; sanity check failed")
          }
        }
      }
    }

    # --------------------------------------------------------------------------
    ## SET TF-CPG PRIOR
    # --------------------------------------------------------------------------
    print("Setting TF-CpG priors.")
    
    # get all TFs
    tfs <- ranges$tfs$SYMBOL
    # also get the chipseq context for our cpgs
    context <- get_tfbs_context(names(ranges$cpgs), fcpgcontext)

    # for all cpgs
    for(c in rownames(context)){
      for(tf in tfs) {
        # might be that the TF was measured in more than one cell line
        if(any(context[c,grepl(tf, colnames(context))])) {
          if(c %in% colnames(priors) &
             tf %in% colnames(priors)){
            # set arbitrary large prior (chip-evidence for this connection..)
            priors[c,tf] <- priors[tf,c] <- 0.99
          } else {
            warning(paste0("WARNING: TF-CpG link not available:",
                           c, "-", tf))
          }
        }
      }
    }
  }

  # ----------------------------------------------------------------------------
  ## SET SNP PRIOR (if available)
  # ----------------------------------------------------------------------------
  # iterate over each of the snp genes and set prior according to
  # pre-calculated priors
  # NOTE: eqtl_priors are already filtered for our SNP
  if(exists("eqtl_priors") && !is.null(eqtl_priors)) {
    print("Setting SNP priors.")
    
    snp_genes <- ranges$snp_genes$SYMBOL
    for(g in snp_genes) {
      # already filtered for sentinel id, just check gene
      idx <- which(grepl(paste0(paste0(",",g,"$"),"|",
                                paste0(",",g,","),"|",
                                paste0("^",g,","),"|",
                                paste0("^",g,"$")), eqtl_priors$symbol))[1]
      # did we find the gene? then set prior
      if(!is.na(idx)){
        p <- (1-eqtl_priors[idx]$lfdr)
        if(p <= pseudo_prior) {
          p <- pseudo_prior
        }
        if(p>=1) {
          p = 1 - pseudo_prior
        }
        priors[g,id] <- priors[id,g] <- p
      }
    }
  } else if(!exists(eqtl_priors)){
    stop("eQTL priors not loaded.")
  }

  # ----------------------------------------------------------------------------
  ## SET GENE-GENE PRIOR
  # ----------------------------------------------------------------------------
  print("Setting Gene-Gene priors.")
  if(!exists("genegene_priors")) {
    stop("genegene priors not loaded!")
  }

  # here we add the priors based on the string interaction network
  genes <- colnames(priors)[!grepl("^rs|^cg", colnames(priors))]

  # get subset of edges which are in our current graph
  PPI_SUB <- subGraph(intersect(nodes(ppi_db), genes), ppi_db)
  sedges <- edgeL(PPI_SUB)
  snodes <- nodes(PPI_SUB)
  rnames <- rownames(genegene_priors)

  temp <- lapply(names(sedges), function(n) {
    l <- sedges[[n]]$edges
    if (length(l) > 0) {
      sapply(l, function(i) {
        m <- snodes[i]
        # set priors for those validated connections found in the PPI DB
        # check whether we have a specific prior for this gene-gene connection
        # match genes either at the beginning with a trailing "_" or
        # at the end with a beginning "_" (rownames are: gene1_gene2)
        idxs <- which((grepl(paste0("^", n, "_"), rnames) &
                         grepl(paste0("_", m, "$"), rnames)) |
                        (grepl(paste0("^", m, "_"), rnames) &
                           grepl(paste0("_", n, "$"), rnames)))[1]
        if(!is.na(idxs)) {
          p <- pseudo_prior + genegene_priors[idxs, "prior"]
          if(p >= 1) p <- 1 - pseudo_prior
        } else {
          p <- pseudo_prior
        }
        return(list(g1=m, g2=n, prior=p))
      })
    }
  })

  temp <- unlist(temp)
  
  # did we have at least one connection?
  if(!is.null(temp)) {
    temp <- matrix(temp, ncol=3, byrow = T)
    for(i in 1:nrow(temp)) {
      # get genes and prior -> set values in prior matrix
      g1 <- temp[i,1]
      g2 <- temp[i,2]
      p <- as.numeric(temp[i,3])
      
      priors[g1, g2] <- priors[g2, g1] <- p
    }
  }

  # sanity check, matrix should not contain 0s or 1s or values smaller than our
  # peudo prior
  if(any(priors==0) | any(priors==1) | any(priors<pseudo_prior)) {
    stop(paste0("ERROR: Sanity check for priors failed. ",
                "(contain 0s/1s or values smaller than ", pseudo_prior, ")"))
  }

  return(priors)
}

#'------------------------------------------------------------------------------
#' Load eqtl priors into environment
#'
#' @param sentinel The sentinel id. eQTL data will immediately
#' be filtered for this SNP id.
#' @param feqtl_priors The RDS file containing the eqtl priors
#' @param prior_type The type of the eQTL priors to be loaded (eqtlgen or gtex)
#' @parma env The environment in which to load the priors. Default: .GlobalEnv
#'
#' @author Johann Hawe
#' @return nothing
#'------------------------------------------------------------------------------
load_eqtl_priors <- function(sentinel,
                             feqtl_priors,
                             prior_type = c("eqtlgen","gtex"),
                             env=.GlobalEnv) {

  # check whether prior were already loaded
  if(with(env, exists("eqtl_priors"))) {
    # nothing to do
    warning("eQTL priors already loaded, doing nothing.")
    return()
  }

  print("Loading eQTL priors.")
  eqtl_priors <- readRDS(feqtl_priors);

  # check whether to filter eqtl data
  if("gtex" %in% prior_type) {
    if(sentinel %in% eqtl_priors$RS_ID_dbSNP135_original_VCF |
       sentinel %in% eqtl_priors$RS_ID_dbSNP142_CHG37p13) {
      eqtl_priors <- eqtl_priors[(eqtl_priors$RS_ID_dbSNP135_original_VCF==sentinel |
                                eqtl_priors$RS_ID_dbSNP142_CHG37p13==sentinel), ]
    } else {
      warning(paste0("Sentinel ", sentinel, " has no GTEx eQTL."))
      eqtl_priors <- NULL
    }
  } else {
    if(sentinel %in% eqtl_priors$SNP) {
      eqtl_priors <- eqtl_priors[eqtl_priors$SNP==sentinel, ]
    } else {
      warning(paste0("Sentinel ", sentinel, " has no GTEx eQTL!"))
      eqtl_priors <- NULL
    }
  }

  assign("eqtl_priors", eqtl_priors, env)
}

#'------------------------------------------------------------------------------
#' Load gene-gene priors into environment
#'
#' @author Johann Hawe
#' @return nothing
#'------------------------------------------------------------------------------
load_genegene_priors <- function(fgene_priors,
                                 env=.GlobalEnv) {
  # check whether prior were already loaded
  if(with(env, exists("genegene_priors"))){
    # nothing to do
    warning("gene-gene priors already loaded, doing nothing.")
    return()
  }

  print("Loading genegene_priors priors.")
  genegene_priors <- readRDS(fgene_priors);

  assign("genegene_priors", genegene_priors, env);
}

#'------------------------------------------------------------------------------
#' Create priors from the GTEX data
#'
#' @author Johann Hawe
#'
#' @return nothing
#'------------------------------------------------------------------------------
create_priors <- function(feqtl, fsnpInfo, frpkm,
                          fsampleDS, fphenotypeDS, dplots, ppi_db,
			  fout_gene_priors, fout_eqtl_priors) {

  # might be the plotting directory does not exist yet
  dir.create(dplots)

  # create the gene-gene priors
  res <- create.gtex.gene.priors(fsampleDS, fphenotypeDS,
                                 frpkm, dplots, ppi_db)
  saveRDS(res, file=fout_gene_priors)

  gc()

  # create the snp-gene priors
  res <- create.gtex.eqtl.priors(feqtl, fsnpInfo)
  saveRDS(res, fout_eqtl_priors)
}

#'------------------------------------------------------------------------------
#' Creates the gtex gene-gene priors under consideration of the
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
                                    rpkm.file, plot.dir, ppi_db) {

  # load the string information
  PPI.NODES <- nodes(ppi_db)
  PPI.EDGES <- graph::edges(ppi_db)

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

  # NOTE: We are not using covariates anymore during normalization since
  # we rely on the first 10 PCs
  gene.matrix <- t(gene.matrix)
  gene.matrix <- log2(gene.matrix+1)
  gene.matrix <- normalize.expression(gene.matrix)
  gene.matrix <- gene.matrix[,which(colnames(gene.matrix) %in% PPI.NODES)]
  # plot gene expression values
  pdf(paste0(plot.dir, "/expression.values.pdf"))
  hist(gene.matrix, breaks=150, xlab="expression residuals")
  dev.off()

  print("Calculating gene-wise correlations.")

  cnames <- colnames(gene.matrix)

  temp <- lapply(names(PPI.EDGES), function(e) {
      lapply(PPI.EDGES[[e]], function(j) {
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
  pdf(paste0(plot.dir, "/expression_correlation.pdf"))
  hist(gtex.gg.cors$correlation, breaks=100, xlab="correlation", main="gtex gene correlations")
  abline(v=0, col="red")
  dev.off()

  # report the pi0 and the proportion of significant tests for qvalue package
  pi1 <- 1-pi0est(gtex.gg.cors[,"pval"])$pi0
  print(paste0("1-pi0 is: ", pi1))

  # return prior
  return(gtex.gg.cors)
}

#'------------------------------------------------------------------------------
#' Preprocesses the gtex Whole_Blood eqtl data to add a column of local FDR
#' values which we will use (1-lfdr) for our snp-gene priors
#' Output is written to a predetermined gz-file
#'
#' @author Johann Hawe
#' @return nothing
#'
#'------------------------------------------------------------------------------
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
  return(pairs)
}

#'------------------------------------------------------------------------------
#' Preprocesses the eQTLgen Whole_Blood eqtl data to add a column of local FDR
#' values which we will use (1-lfdr) for our snp-gene priors
#'
#' @author Johann Hawe
#'
#' @return nothing
#'
#'------------------------------------------------------------------------------
create_eqtlgen_eqtl_priors <- function(eqtl_file) {

  require(data.table)
  require(fdrtool)

  # format: Pvalue  SNP     SNPChr  SNPPos  Zscore  AssessedAllele  OtherAllele     
  # Gene    GeneSymbol      GeneChr GenePos NrCohortsNrSamples       FDR
  pairs <- fread(eqtl_file,
                 header=T,
                 select=c("Pvalue", "SNP", "Gene", "GeneSymbol"),
                 data.table=T)
  colnames(pairs) <- c("Pvalue", "SNP", "Gene", "symbol")

  # get local fdr for all entries
  pairs[, ("lfdr") := fdrtool(pairs$Pvalue, statistic="pvalue")$lfdr]
  pairs[, ("prior") := 1-pairs$lfdr]

  # output rds file
  return(pairs)
}

#'------------------------------------------------------------------------------
#' Method creates priors for TF~CpG links for application of the GGM
#' Loads Remap/ENCODE TFBS and uses score column of bed files to infer
#' prior values.
#'
#' TODO: a possible extension could be to get the sequences where the TFBS are
#' identified to be bound and check their PWM score for those positions...
#'
#' @author Johann Hawe
#'
#'------------------------------------------------------------------------------
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
