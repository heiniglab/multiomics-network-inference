#!/usr/bin/env Rscript

log <- snakemake@log[[1]]
sink(file=log)

library(GenomicRanges)

# parse the sentinels from the available data files
files <- list.files("data/current/cohorts/", ".*ranges.RData", full.names = F)
sentinels <- unlist(lapply(sapply(files, strsplit, "\\."), "[[", 1))
# get data for both cohorts for all sentinels
cohorts <- c("kora", "lolipop")

# for plotting the expression plots
plot.dir <- "results/current/plots/expression_per_sentinel"
if(!dir.exists(plot.dir)){
  dir.create(plot.dir, recursive = T)
}
data <- lapply(sentinels, function(sentinel){
  d <- lapply(cohorts, function(co){
    # load data into new environment
    env <- new.env()
    dfile <- paste0("data/current/cohorts/",
                    sentinel, ".", toupper(co), ".adjusted.data.RData")
    if(!file.exists(dfile)){
      warning(paste0(co, " data does not exist for sentinel ", sentinel, "."))
    } else {
      load(dfile, envir=env)
      ggm.data <- with(env, processed)
      nodes <- colnames(ggm.data)
      # adjust ranges for only those which were initially found in the data
      load(paste0("data/current/cohorts/",
                  sentinel, ".ranges.RData"))
      ranges$cpg.genes <- ranges$cpg.genes[(ranges$cpg.genes$SYMBOL) %in% nodes]
      ranges$snp.genes <- ranges$snp.genes[(ranges$snp.genes$SYMBOL) %in% nodes]
      ranges$spath <- ranges$spath[ranges$spath$SYMBOL %in% nodes]
      ranges$tfs <- ranges$tfs[ranges$tfs$SYMBOL %in% nodes]
      return(list(id=id, data=ggm.data, ranges=ranges, nodes=nodes))
    }
  })
  names(d) <- cohorts
  return(d)
})
names(data) <- sentinels
save(file=snakemake@output[[1]], data)

sink()
