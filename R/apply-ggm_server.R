# this script calls for each of the available sentinels (parsed from the data)
# the ggm modelling scripts.
library(rmarkdown)
args <- commandArgs(trailingOnly = T)
snp <- args[1]
cores <- args[2]

# make sure output dir exists
dir.create("results/reports", recursive = F, showWarnings = F)
render("1-apply-ggm.Rmd", 
       params=list(SNP=snp, cores=cores), 
       output_file = paste0("results/reports/apply-ggm-", snp, ".html"))
