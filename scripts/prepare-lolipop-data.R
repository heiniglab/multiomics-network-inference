# -------------------------------------------------------------------------------
#' Script to collect and preprocess needed lolipop data  for all subsequent
#' analysis
#'
#' @author Johann Hawe
#'
# -------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
print("Get snakemake params.")
# -------------------------------------------------------------------------------
fdata <- snakemake@input[["lolipop"]]
fdata_update <- snakemake@input[["lolipop_update"]]
fout <- snakemake@output[[1]]

# -------------------------------------------------------------------------------
print("Load and prepare data.")
# -------------------------------------------------------------------------------

prepare_data <- function(data_file) {
  load(data_file)

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

  return(list(expr=expr, meth=meth, geno=geno, covars=covars))
}

# prepare the old and update file individually
old <- prepare_data(fdata)
new <- prepare_data(fdata_update)

# ensure that we have the same individuals and ordering
if(!all(rownames(old$geno) == rownames(new$geno))) {
  stop("Individuals/ordering does not match.")
}

# now merge the individual data frames
expr <- cbind(old$expr, new$expr[,!colnames(new$expr) %in% colnames(old$expr)])
geno <- cbind(old$geno, new$geno[,!colnames(new$geno) %in% colnames(old$geno)])

# covars and meth are the same between updates
meth <- old$meth
covars <- old$covars

# just report some stats
print("Old data: -----------")
d <- old
print("Expr dimension:")
print(dim(d$expr))
print("Geno dimension:")
print(dim(d$geno))
print("Meth dimension:")
print(dim(d$meth))
print("Covars dimension:")
print(dim(d$covars))

print("New data: -----------")
d <- new
print("Expr dimension:")
print(dim(d$expr))
print("Geno dimension:")
print(dim(d$geno))
print("Meth dimension:")
print(dim(d$meth))
print("Covars dimension:")
print(dim(d$covars))

print("Merged data: -----------")
print("Expr dimension:")
print(dim(expr))
print("Geno dimension:")
print(dim(geno))
print("Meth dimension:")
print(dim(meth))
print("Covars dimension:")
print(dim(covars))

# ------------------------------------------------------------------------------
print("Saving data.")
# ------------------------------------------------------------------------------
save(file=fout, expr, meth, geno, covars)
