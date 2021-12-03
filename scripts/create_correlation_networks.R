#' -----------------------------------------------------------------------------
#' Creates a correlation based network for a specific sentinel
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Wed May 19 11:15:10 2021
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
source("scripts/lib.R")
source("scripts/reg_net.R")
source("scripts/reg_net_utils.R")

# ------------------------------------------------------------------------------
print("Get snakemake params.")
# ------------------------------------------------------------------------------

fdata_kora <- snakemake@input$kora
fdata_lolipop <- snakemake@input$lolipop
foutput <- snakemake@output$result

# ------------------------------------------------------------------------------
print("Process data.")
# ------------------------------------------------------------------------------

remove_all_na <- function(data) {
  # remove (rare) all-NA cases. This can happen due to scaling of all-zero entities,
  # which can arise due to a very large number of cis-meQTLs which effects get
  # removed from the CpGs during data preprocessing.
  # NOTE: we could possibly handle this differently? Seems that these particular
  # cpgs are highly influenced by genetic effects?
  use <- apply(data, 2, function(x)
    (sum(is.na(x)) / length(x)) < 1)
  data <- data[, use]
  data
}

data_kora <- remove_all_na(readRDS(fdata_kora))
data_lolipop <- remove_all_na(readRDS(fdata_lolipop))

print("Get KORA correlation graph...")
correlation_fit_kora <-
  reg_net(data_kora, NULL, "correlation")
print(correlation_fit_kora$graph)

print("Done.\nGet LOLIPOP correlation graph...")
correlation_fit_lolipop <-
  reg_net(data_lolipop, NULL, "correlation")
print(correlation_fit_lolipop$graph)
print("Done.")

result <- list(kora = correlation_fit_kora,
               lolipop = correlation_fit_lolipop)

readr::write_rds(result, foutput)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
