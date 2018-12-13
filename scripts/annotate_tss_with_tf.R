# ------------------------------------------------------------------------------
#' Creates a TF network annotated with hg19 TSS
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and scripts.")
# ------------------------------------------------------------------------------
library(graph)
library(rtracklayer)
library(data.table)

source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
fbinding_sites_remap <- snakemake@input[["tfbs_remap"]]
fbinding_sites_encode <- snakemake@input[["tfbs_encode"]]

ftfbs_annot <- snakemake@output[["tfbs_annot"]]

# ------------------------------------------------------------------------------
print("Start processing.")
# ------------------------------------------------------------------------------
ga <- get.gene.annotation()
tss <- promoters(ga, 1000, 1000)
names(tss) <- tss$SYMBOL

# ------------------------------------------------------------------------------
print("Creating TF-TSS annotation.")
# ------------------------------------------------------------------------------

#' Creates an annotation object, mapping TFBS to TSS
#'
#' Loads all available TFBS collected from public sources (Encode, Remap) and
#' overlaps those with the provided TSS.
#' Code adapted from file R/annotate-cpgs.R.
#'
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
annotate_tfbs_to_tss <- function(fbinding_sites_remap,
                                 fbinding_sites_encode,
                                 tss) {
  # get the TFBS regions from remap
  tfbs = import(fbinding_sites_remap)
  ann = t(matrix(unlist(strsplit(values(tfbs)[,"name"], ".", fixed=T)), nrow=3))
  colnames(ann) = c("geo_id", "TF", "condition")
  values(tfbs) = DataFrame(name=values(tfbs)[,"name"],
                           data.frame(ann, stringsAsFactors=F))

  # we write out a table with all conditions and select the blood related ones
  conditions = t(matrix(unlist(strsplit(unique(values(tfbs)[,"name"]), ".",
                                        fixed=T)), nrow=3))
  colnames(conditions) = c("geo_id", "TF", "condition")
  conditions = conditions[order(conditions[,"condition"]),]
  conditions = conditions[,c(1,3)]
  conditions = conditions[!duplicated(paste(conditions[,1], conditions[,2])),]
  conditions = data.frame(conditions, blood.related=F)
  for (term in c("amlpz12_leukemic", "aplpz74_leukemia",
                 "bcell", "bjab", "bl41",
                 "blood", "lcl", "erythroid", "gm",
                 "hbp", "k562", "kasumi",
                 "lymphoblastoid", "mm1s", "p493",
                 "plasma", "sem", "thp1", "u937")) {
    conditions[grep(term, conditions[,2]),"blood.related"] = TRUE
  }

  # select the appropriate blood related TFBS subset
  selected = tfbs[values(tfbs)[,"condition"] %in%
                    conditions[conditions[,"blood.related"],"condition"]]

  # load the encode tfs separately
  encode = as.data.frame(fread(fbinding_sites_encode, header=F))
  encode = with(encode, GRanges(seqnames=V1, ranges=IRanges(V2 + 1, V3),
                                name=paste("ENCODE", V4, tolower(V6), sep="."),
                                geo_id="ENCODE", TF=V4,
                                condition=tolower(V6)))

  # filter blood related cell lines
  encode.lcl = encode[grep("gm", values(encode)[,"condition"])]
  values(encode.lcl)[,"condition"] = "lcl"
  encode.k562 = encode[grep("k562", values(encode)[,"condition"])]
  values(encode.k562)[,"condition"] = "k562"

  # combine remap and encode TFBS
  selected = c(selected, encode.lcl, encode.k562)

  # create an annotation matrix for the TSS
  chip = paste(values(selected)[,"TF"], values(selected)[,"condition"], sep=".")
  chip_exp = unique(chip)

  tfbs_ann = sapply(chip_exp, function(x) overlapsAny(tss,
                                                      selected[chip == x]))
  rownames(tfbs_ann) = names(tss)

  return(tfbs_ann)
}

tfbs_annot <- annotate_tfbs_to_tss(fbinding_sites_remap,
                                   fbinding_sites_encode,
                                   tss)
# ------------------------------------------------------------------------------
print("Saving results.")
# ------------------------------------------------------------------------------
saveRDS(tfbs_annot, file=ftfbs_annot)

# ------------------------------------------------------------------------------
print("Session info:")
# ------------------------------------------------------------------------------
sessionInfo()