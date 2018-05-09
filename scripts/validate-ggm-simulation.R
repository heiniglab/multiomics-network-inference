#'
#' Script to validate the GGMs calculated on simulated
#' data
#' 
#' @author Johann Hawe
#'


sink(file=snakemake@log[[1]])

# load needed libraries
library(BDgraph, lib="/storage/groups/groups_epigenereg/packages/2017/R/3.4")
library(GenomicRanges)

source("scripts/lib.R")

ifile <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]

print(paste0("Processing file: ", ifile, "."))

# load data
load(ifile)

# the result table
tab <- c()

result <- result

# check each individual simulated result
temp <- lapply(names(result), function(n) {
  # get validation data
  r <- result[[n]]
  d <- r$data.sim
  ggm_fit <- r$fits$ggm_fit
  ggm_fit_no_priors <- r$fits$ggm_fit_no_priors
  
  # use the bdgraph internal method to get spec/sens, f1 and MCC
  perf <- t(compare(d, ggm_fit, ggm_fit_no_priors))
  comparisons <- c("true", "ggm_fit", "ggm_fit_no_priors")
  perf <- as.data.frame(perf)
  rownames(perf) <- paste(n, comparisons, sep="_")
  # remember for easy plotting
  perf$rdegree <- r$rdegree
  perf$snp <- r$snp
  perf$comparison <- comparisons
  tab <<- rbind(tab, perf)
})

# write result to output file
write.table(file=ofile, tab, col.names=NA, sep="\t", quote=F, row.names = TRUE)

sink()
