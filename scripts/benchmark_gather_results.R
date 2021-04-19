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
theme_set(theme_cowplot() + background_grid())

benchmark_number_iterations <- snakemake@params$benchmark_number_iterations

# ------------------------------------------------------------------------------
print("Loading data")
# ------------------------------------------------------------------------------
result_files <- snakemake@input
data <- lapply(result_files, function(input_file) {
  as_tibble(readRDS(input_file))
}) %>% bind_rows

# overview plot of results

# define colors
paired <- RColorBrewer::brewer.pal(4, "Paired")
names(paired) <- c("glasso", "glasso (priors)", "bdgraph", "bdgraph (priors)")
unpaired <- RColorBrewer::brewer.pal(7, "Dark2")[c(2,3,7)]
names(unpaired) <- c("irafnet", "genie3", "genenet")
graph_cols <- c(paired, unpaired)

# we need slightly nicer names
data <- mutate(data, expr=gsub("_priors"," (priors)", expr))

gp <- ggplot(data,
             aes(x = reorder(expr, -(time)), y = log10(time), color=expr)) +
  scale_y_log10() +
  geom_boxplot() +
  facet_grid(rows = vars(sample_size), cols = vars(number_of_nodes)) + 
  scale_color_manual(values = graph_cols) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0)) +
  theme(panel.border = element_rect(fill = NA, size = 1, color = "lightgrey")) +
  labs(
    x = "",
    y = "log10(time in seconds)",
    title = paste0(
      "Benchmark results for ",
      benchmark_number_iterations,
      " iterations."
    ), 
    subtitle = "Columns show number of input nodes, rows sample size",
    color = "model"
  )
gp

save_plot(filename = snakemake@output$overview_plot, gp, 
          ncol = 2 , nrow = 2)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
devtools::session_info()
