# -------------------------------------------------------------------------------
#' Script to collect and preprocess needed lolipop data  for all subsequent 
#' analysis 
#'
#' @author Johann Hawe
#'
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Get snakemake params
# -------------------------------------------------------------------------------
fdata <- snakemake@input[["lolipop"]]
fout <- snakemake@output[[1]]

# -------------------------------------------------------------------------------
# Load and prepare data
# -------------------------------------------------------------------------------
load(fdata)

# prepare genotypes
geno <- t(dosage.oX_ia)
geno <- geno[,apply(geno,2,var) != 0, drop=F]

meth <- t(beta)
# sort using the genotype ordering
meth <- meth[rownames(geno),]

expr <- t(exma)
# sort using the genotype ordering
expr <- expr[rownames(geno),]
expr <- expr[,!grepl("NA", colnames(expr)),drop=F]

covars <- phe
rownames(covars) <- covars$Sample.ID
covars <- covars[rownames(geno),,drop=F]

# change some colnames to correpond to the same names used in the KORA data
cnames <- colnames(covars)
colnames(covars)[grepl("Sex",cnames)] <- "sex"
colnames(covars)[grepl("Age",cnames)] <- "age"
colnames(covars)[grepl("RNA_conv_batch",cnames)] <- "batch1"
colnames(covars)[grepl("RNA_extr_batch",cnames)] <- "batch2"
covars[,"batch1"] <- factor(covars[,"batch1"])
covars[,"batch2"] <- factor(covars[,"batch2"])

# ------------------------------------------------------------------------------
# All done, save.
# ------------------------------------------------------------------------------
save(file=fout, expr, meth, geno, covars)
