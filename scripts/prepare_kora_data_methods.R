
# ------------------------------------------------------------------------------
#' Get the genotypes for a certain set of snp ranges (width=1)
#'
#' @param snp_ranges A GRanges of SNPs for which to get genotypes. Needed for
#' KORA cohort only
#' @param dosage_file The Tabix indexed file containing the dosages
#' @param individuals The list of individuals (ids) to be found in the dosages
#' file
#' @param individuals_to_keep The individuals to keep
#'
#' @return Matrix of genotypes with SNPs in the columns and individuslasin the 
#' rows
#'
# ------------------------------------------------------------------------------
get_genotypes <- function(snp_ranges, dosage_file,
                          individuals, individuals_to_keep){

  if(is.null(snp_ranges)) {
    return(NULL)
  }

  # split into junks of 100,000 SNPs a piece
  snp_ranges <- split(snp_ranges, ceiling(seq_along(snp_ranges)/100000))
  l <- length(snp_ranges)
  
  print(paste0("Processing ", l, " splits."))
  
  # prepare tabix file reference
  tf <- open(TabixFile(dosage_file))
  temp <- lapply(1:l, function(r) {
    # report progress
    if(r %% 100 == 0) print(paste0("Processed ", r, " splits."))
    
    # get genotypes
    geno <- scan_snps(snp_ranges[[r]], tf, individuals, individuals_to_keep)
    geno
  })
  geno <- do.call(rbind, temp)

  print(warnings())
  
  if(nrow(geno) < l){
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
