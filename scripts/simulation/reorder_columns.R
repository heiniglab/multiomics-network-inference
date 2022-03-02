# ------------------------------------------------------------------
# Used to combine the results of different simulation validation runs
# to ensure proper column ordering over all results.
# ------------------------------------------------------------------

library(magrittr)

fresults <- snakemake@input$results

result <- lapply(fresults, function(f) {
  res <- readr::read_tsv(f) %>%
    dplyr::select(snp, rdegree, comparison, dplyr::everything())
  
  if (grepl("subset", f)) {
    s <- gsub(".*subset([0-9]+).*", "\\1", f)
    res <- dplyr::mutate(res,
                         subset = s)
  }
  res
  
}) %>% dplyr::bind_rows()

readr::write_tsv(result, snakemake@output$combined)
