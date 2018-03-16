#'
#' Script to split the cosmo-RData file into cis/long-range/trans subsets
#'
#' @author Johann Hawe
#'

# get file paths
f_cosmo <- snakemake@input[[1]]
of_trans <- snakemake@output[["trans"]]
of_cis <- snakemake@output[["cis"]]
of_lr <- snakemake@output[["longrange"]]

cat("Loading cosmo file.\n")
load(f_cosmo)

# threshold for long-range: 1MB distance between cpg and snp
cat("Defining cis/longrange/trans.\n")
cis <- which((cosmo$snp.chr == cosmo$cpg.chr) & (abs(cosmo$snp.pos - cosmo$cpg.pos)<=1e6))
cosmo.cis <- cosmo[cis,]
longrange <- which((cosmo$snp.chr == cosmo$cpg.chr) & (abs(cosmo$snp.pos - cosmo$cpg.pos)>1e6))
cosmo.longrange <- cosmo[longrange,]
trans <- which(cosmo$snp.chr != cosmo$cpg.chr)
cosmo.trans <- cosmo[trans,]
rm(cosmo)

cat("Number of cis-associations: ", length(cis), "\n")
cat("Number of longrange-associations: ", length(longrange), "\n")
cat("Number of trans-associations: ", length(trans), "\n")

cat("Saving new cosmo splits.\n")

cosmo <- cosmo.cis
saveRDS(file=of_cis, cosmo)
cosmo <- cosmo.longrange
saveRDS(file=of_lr, cosmo)
cosmo <- cosmo.trans
saveRDS(file=of_trans, cosmo)
