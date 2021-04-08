#' -----------------------------------------------------------------------------
#' Benchmarks individual inference methods with respect to runtime
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Thu Apr  8 07:43:53 2021
#' -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(microbenchmark)
# needed to simulate data
library(BDgraph)

source("scripts/lib.R")
source("scripts/reg_net.R")
source("scripts/benchmark_methods.R")

threads <- snakemake@threads

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

simulation_number_of_nodes <-
  snakemake@params$simulation_number_of_nodes
simulation_sample_size <- snakemake@params$simulation_sample_size
benchmark_number_iterations <-
  snakemake@params$benchmark_number_iterations

model <- snakemake@wildcards$model

# ------------------------------------------------------------------------------
print("Preparing benchmark input functions.")
# ------------------------------------------------------------------------------
benchmark_input <- list()

is_model_with_prior <-
  ifelse(model %in% c("glasso", "bdgraph"),
         TRUE, FALSE)

# only one where we can use  priors but without 'no prior' version
if (model == "irafnet") {
  benchmark_input[[model]] <- function() {
    irafnet = reg_net(s$data, s$priors, model,
                      threads = threads)
  }
  
} else {
  if (is_model_with_prior) {
    benchmark_input[[paste0(model, "_priors")]] <- function() {
      
      reg_net(s$data, s$priors, model,
              threads = threads)
    }
    benchmark_input[[model]] <- function() {
      # we set 'use_gstart' in case it is bdgraph
      reg_net(s$data, NULL, model, use_gstart = F,
              threads = threads)
    }
  } else {
    benchmark_input[[model]] <- function() {
      reg_net(s$data, NULL, model,
              threads = threads)
    }
  }
}

# ------------------------------------------------------------------------------
print("Performing benchmark.")
# ------------------------------------------------------------------------------

benchmark_results <-   microbenchmark(
  list = benchmark_input,
  times = benchmark_number_iterations,
  unit = "s",
  setup = {
    s <-
      simulate_data(simulation_number_of_nodes,
                    simulation_sample_size)
  }
)

# ------------------------------------------------------------------------------
print("Benchmark done. Finishing up.")
# ------------------------------------------------------------------------------
saveRDS(benchmark_results, file = snakemake@output$result_table)

# ------------------------------------------------------------------------------
print("Report warnings:")
# ------------------------------------------------------------------------------
warnings()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
devtools::session_info()
