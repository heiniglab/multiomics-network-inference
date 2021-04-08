#' -----------------------------------------------------------------------------
#' Gathers all results from the benchmark analysis and creates an overview plot
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Apr  8 11:25:30 2021
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(microbenchmark)
library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot() + background_grid(major = "xy"))

# ------------------------------------------------------------------------------
print("Loading data")
# ------------------------------------------------------------------------------
data <- lapply(snakemake@input, function(input_file) {
  as_tibble(readRDS(input_file))
}) %>% bind_rows

# quick overview plot of results
gp <- ggplot(data,
             aes(x = reorder(expr, -(time)), y = log10(time))) +
  scale_y_log10() +
  geom_boxplot() +
  labs(
    x = "model",
    y = "log10(time in seconds)",
    title = paste0(
      "Benchmark results for ",
      benchmark_number_iterations,
      " iterations."
    ),
    subtitle = paste0(
      "P=",
      simulation_number_of_nodes,
      ";  N=",
      simulation_sample_size,
      ";  threads=",
      threads
    )
  )

save_plot(filename = snakemake@output$overview_plot, gp)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
devtools::session_info()