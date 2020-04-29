library(GenomicRanges)
library(GenomicFeatures)
library(FDb.InfiniumMethylation.hg19)
library(data.table)
library(illuminaHumanv3.db)
library(rtracklayer)
library(graph)
library(RBGL) # for shortest paths
library(Matrix)
library(xtable)

source("scripts/lib.R")
source("scripts/collect_ranges_methods.R")

missing <- read_tsv("missing_sentintels.txt")$sentinel

# for the ranges collection
fcosmo <- "results/current/trans-cosmopairs_combined_151216.rds"
fmeqtl <- "data/current/meQTLs/transpairs_r02_110117_converted_1MB.txt"
fppi_db <- "results/current/biogrid.rds"
fprio_tab <- "data/current/rw_string_v9_ld_wb_prioritize_full_with_empirical_p_lte_0.05_bck.txt"
fgene_annot <- "data/current/gencode_annotations/gencode.v19.annotation.gene.gtf"
fcpgcontext <- "results/current/cpg_context.rds"

gene_annot <- load_gene_annotation(fgene_annot)
gene_annot$ids <- probes.from.symbols(gene_annot$SYMBOL,
                                      as_list=T)
ppi_db <- readRDS(fppi_db)

# load trans-meQTL table
trans_meQTL = read.csv(fmeqtl,
                       sep="\t", stringsAsFactors=F)

# load trans-cosmo information
cosmo <- readRDS(fcosmo)

# load priorization table
prio <- read.table(fprio_tab, sep="\t", header=T, stringsAsFactors = F)


# ------------------------------------------------------------------------------
# Load and preprocess data
print("Loading data.")

info <- lapply(missing, function(sentinel) {
  
  print(paste0("Sentinel is ", sentinel))
  
  # get trans-genes which should be used for shortest path extraction
  sprio <- prio[prio$sentinel == sentinel,,drop=F]
  if(nrow(sprio) > 0) {
    best_trans <- unique(sprio$trans.gene)
  } else {
    best_trans <- NULL
  }
  
  pairs = which(trans_meQTL[,"sentinel.snp"] == sentinel)
  
  # get sentinel idx
  sentinel_idx <- which(cosmo$snp==sentinel)[1]
  
  # the large interval range for the sentinel
  chr <- paste0("chr", trans_meQTL[pairs,"chr.snp"][1])
  start <- trans_meQTL[pairs,"interval.start.snp"][1]
  end <- trans_meQTL[pairs,"interval.end.snp"][1]
  
  sentinel_extrange <- GRanges(chr, IRanges(start,end))
  sentinel_range <- with(cosmo[sentinel_idx,],
                         GRanges(paste0("chr", snp.chr),
                                 IRanges(snp.pos, width=1)))
  names(sentinel_range) <- names(sentinel_extrange) <- sentinel
  
  # get cosmo subset
  idxs <- get.trans.cpgs(sentinel, trans_meQTL, cosmo)
  
  # get related genes, i.e. genes near meQTL loci (snp+cpg)
  
  # extend cpgs
  cosmosub <- cosmo[idxs,]
  
  croi <- with(cosmosub, GRanges(paste0("chr", cpg.chr),
                                 IRanges(cpg.pos,width=2)))
  names(croi) <- as.character(cosmosub[,"cpg"])
  croi <- unique(croi)
  
  # extended sentinel region
  sroi <- sentinel_extrange
  names(sroi) <- sentinel
  
  # get the relevant snp genes by overlapping with our sentinel region
  genes_sroi <- subsetByOverlaps(gene_annot, sroi, ignore.strand=T)
  
  if(length(unlist(genes_sroi$ids)) == 0) {
    print("No IDs for SNP genes.")
    return(as.data.frame(genes_sroi))
  }
})

# generate table of missing sentinels for manuscript
names(info) <- missing
df <- bind_rows(info, .id="sentinel") %>% 
  as_tibble() %>%
  dplyr::select(sentinel, chr=seqnames, cis_gene=SYMBOL, gene_start=start,
         gene_end = end, gene_strand=strand, gene_biotype=BIOTYPE)

df %>% xtable(caption = "Sentinels and their annotated cis genes removed from analysis due to the genes not being measured on the microarrays. Sentinels rs1570038 and rs7924137 did not have any cis genes annotated.", 
              label = "stab:filtered_sentinels")
