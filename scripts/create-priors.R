#' Script to create the GTEX priros (gene-gene and snp-gene priors)
#'
#' @author Johann Hawe

library(qvalue)
library(data.table)
library(graph)
library(fdrtool)
library(Homo.sapiens)

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

print("Loading string db.")
string_db <- readRDS(fstring)

# simply delegate
create.priors(feqtl, fsnpinfo, frpkm, fsampleDS, fphenotypeDS, dplots, string_db)
