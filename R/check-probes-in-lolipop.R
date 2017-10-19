# collect all currently assembled expression probes

files <- list.files("~/work/analysis/meQTLs/results/current/ggm/sentinel_only/correlations/KORA/shortest_paths/", "*.collected.ranges.RData", full.names =T)

exprobes <- lapply(files, function(f){
   load(f)
   cr <- collected.ranges
   return(c(cr$snp.genes$ids, cr$cpg.genes$ids, cr$tfs$ids, cr$spath$ids))  
})

exprobes <- unique(unlist(exprobes))

load("~/epigenereg/datasets/kora/meQTLs/20160816/ggmdata_200217.RData")
cat(paste0(exprobes[!exprobes %in% rownames(exma)],collapse = "\n"))