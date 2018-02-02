#' Script to create the GTEX priros (gene-gene and snp-gene priors)
#'
#' @author Johann Hawe

# source necessary scripts
source("R/lib.R")
source("R/priors.R")

# TODO get filenames etc from snakemake variable
# simply delegate...
create.priors()
