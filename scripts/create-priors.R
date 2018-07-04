#' Script to create the GTEX priros (gene-gene and snp-gene priors)
#'
#' @author Johann Hawe

library(qvalue)
library(data.table)
library(graph)
library(parallel)
library(fdrtool)
library(Homo.sapiens)
library(rtracklayer)
library(FDb.InfiniumMethylation.hg19)

# source necessary scripts
source("scripts/lib.R")
source("scripts/priors.R")

feqtl <- snakemake@input[["eqtl"]]
fsnpinfo <- snakemake@input[["snpinfo"]]
fexpr <- snakemake@input[["expr"]]
fsampleinfo <- snakemake@input[["sampleinfo"]]
fpheno <- snakemake@input[["pheno"]]
fstring <- snakemake@input[["string"]]
dplots <- snakemake@params$plot_dir

threads <- snakemake@threads

print("Loading string db.")
string_db <- readRDS(fstring)

# simply delegate
create.priors(feqtl, fsnpinfo, frpkm, fsampleDS, fphenotypeDS, dplots, string_db)

# ------------------------------------------------------------------------------
# prepare the Banovich based priors, i.e. TF-CpG priors
# ------------------------------------------------------------------------------
# methylation data
meth <- fread("data/current/banovich-2017/methylation/full_matrix.txt", data.table=F)
rownames(meth) <- meth$V1
meth$V1 <- NULL
cpgs <- features(FDb.InfiniumMethylation.hg19)
cpgs <- cpgs[rownames(meth)]

# expression data
expr <- read.table("data/current/banovich-2017/xun_lan/allTFexp.withHeader", 
                   header=T, 
                   sep="\t",
                   stringsAsFactors=F)

# apparently the table contains duplicated entries, remove them
expr <- expr[!duplicated(expr),]
rownames(expr) <- unique(expr[,1])

samples <- intersect(colnames(expr), colnames(meth))

expr <- t(expr[,samples])
meth <- t(meth[,samples])

# get (our) chip-seq context for the cpgs
tfbs_ann <- get.chipseq.context(names(cpgs), fcpgcontext)

pairs <- lapply(colnames(expr), function(tf) {
  # get columns for tf 
  sub <- tfbs_ann[,grepl(tf, colnames(tfbs_ann), ignore.case = T), drop=F]
  rs <- rowSums(sub)
  bound_cpgs <- names(rs[rs>0])
  
  assoc <- unlist(mclapply(bound_cpgs, function(c) {
                         cor.test(expr[,tf], 
                             meth[,c], 
                             method="pearson")$p.value
    }, mc.cores=threads))
  
  cbind.data.frame(TF=rep(tf, length(assoc)), 
                   CpG=bound_cpgs, 
                   rho=assoc, 
                   stringsAsFactors=F)
})

tab <- do.call(rbind, pairs)
colnames(tab) <- c("TF", "CpG", "pval")
tab$qval <- qvalue(tab$pval)$lfdr
tab$prior <- 1-tab$qval
head(tab)
write.table(file="results/current/tf-cpg-prior.txt", sep="\t", col.names=NA,
            row.names=T,
            quote=F, tab)

