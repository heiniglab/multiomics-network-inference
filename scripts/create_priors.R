# ------------------------------------------------------------------------------
#' Script to create the GTEX priors (gene-gene and snp-gene priors)
#'
#' @author Johann Hawe
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# ------------------------------------------------------------------------------
# Load libraries and source scripts
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(graph))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(fdrtool))
suppressPackageStartupMessages(library(Homo.sapiens))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(FDb.InfiniumMethylation.hg19))
source("scripts/lib.R")
source("scripts/priors.R")

# ------------------------------------------------------------------------------
# Get snakemake params
# ------------------------------------------------------------------------------

# inputs
feqtl <- snakemake@input[["eqtl"]]
fsnpinfo <- snakemake@input[["snpinfo"]]
fexpr <- snakemake@input[["expr"]]
fsampleinfo <- snakemake@input[["sampleinfo"]]
fpheno <- snakemake@input[["pheno"]]
fppi <- snakemake@input[["ppi"]]
dplots <- snakemake@params$plot_dir

# outputs
fout_gene_priors <- snakemake@output$gene_priors
fout_eqtl_priors <- snakemake@output$eqtl_priors

# ------------------------------------------------------------------------------
# Start processing
# ------------------------------------------------------------------------------

print("Loading PPI db.")
ppi_db <- readRDS(fppi)

# simply delegate
create_priors(
  feqtl,
  fsnpinfo,
  fexpr,
  fsampleinfo,
  fpheno,
  dplots,
  ppi_db,
  fout_gene_priors,
  fout_eqtl_priors
)

if (FALSE) {
  # ------------------------------------------------------------------------------
  print("Prepare the Banovich based priors, i.e. TF-CpG priors.")
  # ------------------------------------------------------------------------------
  # methylation data
  meth <-
    fread("data/current/banovich-2017/methylation/full_matrix.txt",
          data.table = F)
  rownames(meth) <- meth$V1
  meth$V1 <- NULL
  cpgs <- features(FDb.InfiniumMethylation.hg19)
  cpgs <- cpgs[rownames(meth)]

  # expression data
  expr <-
    read.table(
      "data/current/banovich-2017/xun_lan/allTFexp.withHeader",
      header = T,
      sep = "\t",
      stringsAsFactors = F
    )

  # apparently the table contains duplicated entries, remove them
  expr <- expr[!duplicated(expr), ]
  rownames(expr) <- unique(expr[, 1])

  samples <- intersect(colnames(expr), colnames(meth))

  expr <- t(expr[, samples])
  meth <- t(meth[, samples])

  # ------------------------------------------------------------------------------
  print("Get (our) chip-seq context for the cpgs.")
  # ------------------------------------------------------------------------------
  tfbs_ann <- get_tfbs_context(names(cpgs), fcpgcontext)

  # ------------------------------------------------------------------------------
  print("For each TF, get the correlation to each of the CpGs it is bound nearby")
  # ------------------------------------------------------------------------------
  pairs <- lapply(colnames(expr), function(tf) {
    # get columns for tf
    sub <-
      tfbs_ann[, grepl(tf, colnames(tfbs_ann), ignore.case = T), drop = F]
    rs <- rowSums(sub)
    bound_cpgs <- names(rs[rs > 0])

    assoc <- unlist(mclapply(bound_cpgs, function(c) {
      cor.test(expr[, tf],
               meth[, c],
               method = "pearson")$p.value
    }, mc.cores = threads))

    cbind.data.frame(
      TF = rep(tf, length(assoc)),
      CpG = bound_cpgs,
      rho = assoc,
      stringsAsFactors = F
    )
  })

  # ------------------------------------------------------------------------------
  print("Collect and finalize results.")
  # ------------------------------------------------------------------------------
  tab <- do.call(rbind, pairs)
  colnames(tab) <- c("TF", "CpG", "pval")
  tab$qval <- qvalue(tab$pval)$lfdr
  tab$prior <- 1 - tab$qval
  head(tab)
  write.table(
    file = "results/current/tf-cpg-prior.txt",
    sep = "\t",
    col.names = NA,
    row.names = T,
    quote = F,
    tab
  )
}
# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
