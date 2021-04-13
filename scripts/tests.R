#' -----------------------------------------------------------------------------
#' Contains some quick tests to check data and sanity assumptions
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Tue Apr 13 09:02:34 2021
#' -----------------------------------------------------------------------------

test_node_difference_between_cohorts <- function() {
  # list all data files
  files_lolipop <-
    list.files(
      "results/current/biogrid_stringent/cohort_data_expr/lolipop/",
      ".*_meqtl.rds",
      full.names = T
    )
  files_kora <-
    list.files(
      "results/current/biogrid_stringent/cohort_data_expr/kora/",
      ".*_meqtl.rds",
      full.names = T
    )
  # get rid of raw data (contains also PCs etc.)
  files_lolipop <- files_lolipop[!grepl("raw", files_lolipop)]
  files_kora <- files_kora[!grepl("raw", files_kora)]
  
  # get and name by sentinels
  sentinels <- gsub("_.*", "", basename(files_lolipop))
  names(files_lolipop) <- sentinels
  sentinels <- gsub("_.*", "", basename(files_kora))
  names(files_kora) <- sentinels
  
  # get differences between number of available variables/nodes
  differences <-
    unlist(lapply(sentinels, function(s) {
      a <-
        readRDS(files_kora[[s]])
      b <- readRDS(files_lolipop[[s]])
      abs(ncol(a) - ncol(b))
    }))
  summary(differences)
  table(differences)
}
