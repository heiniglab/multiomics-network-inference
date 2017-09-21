# this script calls for each of the available sentinels (parsed from the data)
# the ggm modelling scripts.
library(rmarkdown)
snp <- commandArgs(trailingOnly = T)

#snp <- "rs9859077"
snp <- "rs7783715"

# make sure output dir exists
dir.create("results/reports", recursive = F, showWarnings = F)
render("1-apply-ggm.Rmd", 
       params=list(SNP=snp), 
       output_file = paste0("results/reports/apply-ggm-", snp, ".html"))
