# ------------------------------------------------------------------------------
#' Script to collect needed kora data for all subsequent analysis in a single
#' RData file. Makes collecting data for sentinel SNPs easier.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts.")
# ------------------------------------------------------------------------------
library(GenomicRanges)
library(Rsamtools)
library(data.table)
source("scripts/lib.R")
source("scripts/biomaRt.R")

# ------------------------------------------------------------------------------
# Method definitions
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#' Method to get the genotypes for a certain set of SNPs
#'
#' @param snp_ranges A GRanges of SNPs for which to get genotypes. Needed for
#' KORA cohort only
#'
#' @return Matrix of genotypes with SNPs in the columns and subjects in the rows
#'
# ------------------------------------------------------------------------------
get_genotypes <- function(snp_ranges, dosage_file,
                          individuals, individuals_to_keep, threads){

  if(is.null(snp_ranges)) {
    return(NULL)
  }

  # split into junks
  snp_ranges <- split(snp_ranges, ceiling(seq_along(snp_ranges)/10000))

  # prepare tabix file reference
  tf <- open(TabixFile(dosage_file))
  temp <- mclapply(snp_ranges, function(r) {
    geno <- scan_snps(r, tf, individuals)
    if(!is.null(geno)) {
      geno <- geno[, individuals_to_keep,drop=F]
    }
    geno
  }, mc.cores = threads)
  geno <- do.call(rbind, temp)

  if(nrow(geno) < length(snp_ranges)){
    warning("Some SNPs were NA in genotype data.")
  }

  geno <- t(geno)

  # get rid of non-changing snps
  if(any(apply(geno,2,var)==0)) {
    warning("Removing non-varying SNPs.")
    geno <- geno[,apply(geno,2,var)!=0, drop=F]
  }
  
  # refactor column names (get rid of beginning "1.")
  colnames(geno) <- gsub("^[0-9]+\\.", "", colnames(geno))
}

# start processing -------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

# inputs
fdosage <- snakemake@input[["genotypes"]]
findividual_mapping <- snakemake@input[["individuals"]]
fexpression <- snakemake@input[["expression"]]
fexpression_cov <- snakemake@input[["expression_cov"]]
fmethylation <- snakemake@input[["methylation"]]
fmethylation_cov <- snakemake@input[["methylation_cov"]]
ftrans_meqtl <- snakemake@input[["trans_meqtl"]]
fhouseman <- snakemake@input[["houseman"]]
fcosmo <- snakemake@input[["cosmo"]]
fceqtl <- snakemake@input[["kora_ceqtl"]]
feqtlgen <- snakemake@input[["eqtl_gen"]]
fimpute_indiv <- snakemake@input[["impute_indiv"]]

# params
threads <- snakemake@threads

# ------------------------------------------------------------------------------
print("Load the imputation individuals (individuals with genotypes).")
# ------------------------------------------------------------------------------
imputation_individuals <- read.table(fimpute_indiv,
                                     header=F)[,1]

# ------------------------------------------------------------------------------
print("Loading snp ranges.")
# ------------------------------------------------------------------------------
trans_meqtl <- read.table(ftrans_meqtl, sep="\t", header=T)
ranges <- unique(with(trans_meqtl,
                      GRanges(paste0("chr", chr.snp),
                              IRanges(pos.snp, width=1))))

ceqtl <- read.table(fceqtl, header=T, sep=";")
ceqtl <- get_snpPos_biomart(ceqtl[,1])
ceqtl <- ceqtl[!grepl("^HG|^HS", ceqtl$chr),]
ranges <- unique(c(ranges, with(ceqtl,
                                GRanges(paste0("chr", chr),
                                        IRanges(start, width=1)))))

eqtlgen <- fread(paste0("zcat ", feqtlgen))
eqtlgen <- unique(with(eqtlgen, GRanges(paste0("chr", SNPChr),
                                        IRanges(SNPPos, width=1))))
ranges <- unique(c(ranges, eqtlgen))

load(fcosmo)
ranges <- unique(c(ranges, with(cosmo,
                                GRanges(paste0("chr", snp.chr),
                                        IRanges(snp.pos, width=1)))))
rm(cosmo, ceqtl)
gc()

ranges <- sort(ranges)
print(paste0("Got ", length(ranges), " genotype ranges to process."))

# ------------------------------------------------------------------------------
print("Loading KORA data.")
# ------------------------------------------------------------------------------

# load id mapping to intersect expr/meth/geno data
id_map <- read.table(findividual_mapping,
                     header=T, sep=";", stringsAsFactors = F)

# convert ids to character (instead of int) for easier
# subsetting of data matrices
id_map$axiom_s4f4 <- as.character(id_map$axiom_s4f4)
id_map$genexp_s4f4ogtt <- as.character(id_map$genexp_s4f4ogtt)
id_map$methylierung_f4 <- as.character(id_map$methylierung_f4)
# drop individuals with NAs (e.g. due to missing BMI)
id_map <- na.omit(id_map)

# ------------------------------------------------------------------------------
print("Preparing KORA raw data.")
# ------------------------------------------------------------------------------

# load covariates to get ids of samples with expression/methylation data
load(fexpression_cov)
expr_sids <-as.character(covars.f4$ZZ.NR)
load(fmethylation_cov)
meth_sids <- rownames(pcs)

# gets us 687 individuals, having all data available (some ids of the id
# map are not contained within the data frame...)
toUse <- which(id_map$genexp_s4f4ogtt %in% expr_sids &
                 id_map$methylierung_f4 %in% meth_sids &
                 id_map$axiom_s4f4 %in% imputation_individuals)

id_map <- id_map[toUse,]
id_map$utbmi <- NULL
id_map$ul_wbc <- NULL
id_map$utalteru <- NULL

print(paste0("Using ", nrow(id_map), " samples."))

# ------------------------------------------------------------------------------
print("Load genotypes.")
# ------------------------------------------------------------------------------
geno <- get_genotypes(ranges, fdosage,
                      imputation_individuals, id_map$axiom_s4f4, threads)
print(paste0("Genotype dimensions: ", paste(dim(geno), collapse=",")))
gc()

# ------------------------------------------------------------------------------
print("Load and process expression and methylation data incl. covariates.")
# ------------------------------------------------------------------------------
load(fexpression)
load(fmethylation)

# sort our input data s.t. each row corresponds to the same individual
# using the created ID mapping table
meth <- t(beta[,id_map$methylierung_f4])
# get rid of zero-variance probes
meth <- meth[,apply(meth,2,var, na.rm=T)!=0, drop=F]
print(paste0("Meth dimensions: ", paste(dim(meth), collapse=",")))
# get the methylation PCA results
pcs <- pcs[id_map$methylierung_f4,]
colnames(pcs) <- paste(colnames(pcs), "cp", sep="_")
meth.pcs <- pcs
print(paste0("Meth cov dimensions: ", paste(dim(meth.pcs), collapse=",")))
rm(pcs)

# load houseman blood count data for methylation
houseman <- read.table(fhouseman,
                       sep=";", header=T,row.names=1)
houseman <- houseman[id_map$methylierung_f4,]
print(paste0("Houseman dimensions: ", paste(dim(houseman), collapse=",")))

# use only those individuals for which we have all data available
expr <- t(f4.norm[,id_map$genexp_s4f4ogtt])
# remove zero-variance probes...
expr <- expr[,apply(expr,2,var)!=0, drop=F]
print(paste0("Expr dimensions: ", paste(dim(expr), collapse=",")))

# load technical covariates for expression data
rownames(covars.f4) <- covars.f4[,1]
covars.f4 <- covars.f4[id_map$genexp_s4f4ogtt,c(2:6)]
covars.f4$sex <- as.factor(covars.f4$sex)
print(paste0("Expression cov dimensions: ",
             paste(dim(covars.f4), collapse=",")))

# create initial data frame with covariates
covars <- cbind(as.data.frame(covars.f4), houseman, meth.pcs)
rm(covars.f4)
gc()

# replace some colnames to easily match with the lolipop data
cc <- colnames(covars)
cc[grepl("storage.time",cc)] <- "batch2"
cc[grepl("plate",cc)] <- "batch1"
colnames(covars) <- cc
covars[,"batch1"] <- factor(covars[,"batch1"])

# ------------------------------------------------------------------------------
print("Saving data.")
# ------------------------------------------------------------------------------
ofile <- snakemake@output[[1]]
save(file=ofile, expr, meth, geno, covars)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()
