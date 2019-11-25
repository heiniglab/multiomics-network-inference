#' -----------------------------------------------------------------------------
#' Compare old (expression based) and new (tfa based) GGM validation results
#'
#' @author Johann Hawe <johann.hawe@helmholtz-muenchen.de>
#'
#' @date Mon Oct 21 18:22:43 2019
#' -----------------------------------------------------------------------------

# ------------------------------------------------------------------------------
print("Load libraries and source scripts")
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

#' Plots several plots comparing TFA and expression
#' results
plot_overview <- function(data) {
  # perform actual plotting
  theme_update(axis.title.x=element_blank(), axis.text.x=element_text(size=12, angle=45, hjust = 1))
  
  a <- data %>%
    ggplot(aes(x=graph_type, y=number_edges, color=type)) + 
    geom_boxplot(position = "dodge") +
    labs(title="number edges")
  
  b <- data %>%
    ggplot(aes(x=graph_type, y=graph_density, color=type)) + 
    geom_boxplot(position = "dodge") + 
    labs(title="graph_density")
  
  c <- data %>%
    ggplot(aes(x=graph_type, y=graph_score, color=type)) + 
    geom_boxplot(position = "dodge") + 
    labs(title="graph_score")
  
  d <- data %>%
    ggplot(aes(x=graph_type, y=cross_cohort_mcc, color=type)) + 
    geom_boxplot(position = "dodge") + 
    labs(title="cross cohort mcc")
  
  e <- data %>%
    ggplot(aes(x=graph_type, y=cross_cohort_f1, color=type)) + 
    geom_boxplot(position = "dodge") +
    labs(title="cross cohort F1 score")
  
  f <- data %>%
    ggplot(aes(x=graph_type, y=cluster, color=type)) + 
    geom_boxplot(position = "dodge") + 
    labs(title="cluster")
  
  cowplot::plot_grid(a, b, c, d, e, f,
                     ncol=2, labels = "AUTO")
}

# read the validation results for meqtls
expr <- read_tsv("results/current/biogrid/validation_expr/validation_all_meqtl.txt")
tfa <- read_tsv("results/current/biogrid/validation_tfa/validation_all_meqtl.txt")

expr <- expr %>% filter(sentinel %in% tfa$sentinel)

data <- bind_rows(expr,tfa)
data <- mutate(data, type=c(rep("expr", nrow(expr)), rep("tfa", nrow(tfa))))

plot_overview(data)

# do the same for the eQTLgen results
expr <- read_tsv("results/current/biogrid_stringent/validation_expr/validation_all_eqtlgen.txt")
tfa <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_eqtlgen.txt")

expr <- expr %>% filter(sentinel %in% tfa$sentinel)

data <- bind_rows(expr,tfa)
data <- mutate(data, type=c(rep("expr", nrow(expr)), rep("tfa", nrow(tfa))))

plot_overview(data)

# ------------------------------------------------------------------------------
print("Done.\nSessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
