#' -----------------------------------------------------------------------------
#' Create a graph object reflecting BIOGRID PPIs which can be used for
#' further analysis.
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Nov 12 16:22:00 2018
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(data.table)
library(graph)
library(igraph)
library(dplyr)
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------
fstring <- snakemake@input$stringdb
fbiogrid <- snakemake@input$biogrid
fhuri <- snakemake@input$huri

fgene_annot <- snakemake@input$gene_annot

# gtex gene expression data (we filter genes for the ones expressed in blood)
fgtex <- snakemake@input$gtex

fout <- snakemake@output[[1]]

ppi_name <- snakemake@params$ppi_name

# ------------------------------------------------------------------------------
print(paste0("Loading PPI database: ", ppi_name))
# ------------------------------------------------------------------------------

ga <- load_gene_annotation(fgene_annot)

ga_tibble <- ga %>% as.data.frame %>% as_tibble(rownames = "ENSG") %>%
  mutate(ENSG_gene = gsub("\\..*", "", ENSG)) %>%
  select(ENSG_gene, SYMBOL)

if(ppi_name == "string") {
  string.all <- fread(fstring,
                      data.table=F, header=T, stringsAsFactors=F)
  string.inter <- string.all[string.all$experimental>=1 | string.all$database>=1,]

  string.nodes <- unique(c(string.inter[,1], string.inter[,2]))
  print("Creating interaction graph.")
  ppi_db <- graphNEL(nodes=string.nodes)
  ppi_db <- addEdge(string.inter[,1],
                    string.inter[,2],
                    ppi_db)
} else if (grepl("biogrid", ppi_name)) {
  # biogrid database
  biogrid <- fread(fbiogrid, stringsAsFactors=F)
  if(grepl("stringent", ppi_name)) {
    print("Making PPIs more stringent.")
    biogrid <- biogrid %>%
      filter(grepl("Low Throughput", Throughput)) %>%
      filter(`Experimental System Type` == "physical") %>%
      as_tibble()
  } else {
    biogrid <- biogrid %>% as_tibble()
  }

  print("Creating interaction graph.")
  nodes <- unique(c(biogrid$`Official Symbol Interactor A`,
                  biogrid$`Official Symbol Interactor B`))
  ppi_db <- graphNEL(nodes=nodes)
  ppi_db <- addEdge(biogrid$`Official Symbol Interactor A`,
                    biogrid$`Official Symbol Interactor B`,
                    ppi_db)
} else if (ppi_name == "huri"){
  huri <- readr::read_tsv(fhuri, col_names = c("gene1", "gene2")) %>% 
    left_join(ga_tibble, by=c("gene1" = "ENSG_gene")) %>% 
    left_join(ga_tibble, by = c("gene2" = "ENSG_gene")) %>% 
    select(gene1 = SYMBOL.x, gene2 = SYMBOL.y) %>% 
    tidyr::drop_na()
  nodes <- unique(c(pull(huri, gene1), pull(huri, gene2)))
  ppi_db <- graphNEL(nodes = nodes)
  ppi_db <- addEdge(huri$gene1, huri$gene2, ppi_db)
} else {
  stop(paste0("PPI db not supported:", ppi_name))
}

# ------------------------------------------------------------------------------
print("Filtering PPI for expressed genes.")
# ------------------------------------------------------------------------------
expr <- fread(fgtex, stringsAsFactors=F)
expressed <- unlist(expr[`Whole Blood` > 0.1,"Name"])
expressed.symbols <- ga[intersect(expressed, names(ga))]$SYMBOL

nodes_to_keep = intersect(nodes(ppi_db), c(expressed.symbols))
ppi_db_expr = subGraph(nodes_to_keep, ppi_db)

# ------------------------------------------------------------------------------
print("Getting largest connected component.")
# ------------------------------------------------------------------------------
ig = graph_from_graphnel(ppi_db_expr)
cl = clusters(ig)
keep = nodes(ppi_db_expr)[cl$membership == which.max(cl$csize)]
ppi_db_final = subGraph(keep, ppi_db_expr)


print("Largest CC in PPI_DB filtered for blood expressed genes:")
print(ppi_db_final)

# ------------------------------------------------------------------------------
print("All done. Saving output.")
# ------------------------------------------------------------------------------
saveRDS(ppi_db_final, fout)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
