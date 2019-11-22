#' -----------------------------------------------------------------------------
#' Script to generate the individual plots we use in the manuscript.
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
library(RColorBrewer)
library(ggpubr)
library(scales)

# set up theme and colors
theme_set(theme_cowplot())
theme_update(legend.text = element_text(size=10), 
             legend.title=element_text(size=10))

sfb_graphs <- scale_fill_brewer(palette="Set2")
scb_graphs <- scale_color_brewer(palette="Set2")
sfb_binary <- scale_fill_brewer(palette = "Paired")
scb_binary <- scale_color_brewer(palette = "Paired")
bgm <- background_grid(major = "xy")
group_cols <- brewer.pal("Set2", n=3)
COLORS <- list(MEQTL = group_cols[1],
               EQTL = group_cols[2])

# ------------------------------------------------------------------------------
print("Figure 1 - Panel A")
# ------------------------------------------------------------------------------
# input directory containing individual simulation validation
# for the meQTL analysis
dinput <- "results/current/biogrid_stringent/simulation/validation/"
finput <- list.files(dinput, "*.txt", full.names = T)

# create data-matrix
temp <- lapply(finput, function(f) {
  f <- read_tsv(f)
  if(nrow(f) > 0) {
    f <- f %>%
      mutate(R = paste0("R=", rdegree)) %>%
      rename(name = X1)
  }
})
tab <- bind_rows(temp)

# create nicer method names
tab <- tab %>% 
  mutate(comparison = gsub("bdgraph$", "bdgraph (priors)", comparison),
         comparison = gsub("glasso$", "glasso (priors)", comparison),
         comparison = gsub("_no_priors", "", comparison))

# get the MCC plot
simulation_mcc <- ggplot(tab,
                         aes(y=MCC, 
                             x=R, 
                             color=reorder(comparison, -MCC, median))) +
  stat_boxplot(geom="errorbar", width=.75)+
  geom_boxplot(outlier.size=0, alpha=0.5, coef=0, outlier.shape = NA) + 
  stat_summary(fun.y=median, geom="smooth", 
               position=position_dodge(0.75),
               aes(group=comparison),lwd=0.8) +
  scb_graphs +
  geom_boxplot(data = tab, alpha=0.5, aes(y=density_true, x=R), 
               fill="#666666", color="#666666",
               inherit.aes = FALSE,
               width=.3) +
  scale_y_continuous(limits=c(min(tab$MCC),1), 
                     sec.axis = sec_axis(trans = ~ ., 
                                         name="true graph density",
                                         breaks=seq(0,0.7,by=0.1))) +
  background_grid(major="xy") +
  labs(x="prior error",
       y="MCC",
       fill="", color="method") + 
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
print("Figure 1 - Panel C")
# ------------------------------------------------------------------------------
# read the validation results for meqtls and eqtls and combine them
meqtl_expr <- read_tsv("results/current/biogrid_stringent/validation_expr/validation_all_meqtl.txt")
meqtl_tfa <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_meqtl.txt")
meqtl <- bind_rows(meqtl_expr,meqtl_tfa) %>%
  mutate(type=c(rep("expr", nrow(meqtl_expr)), rep("tfa", nrow(meqtl_tfa))),
         qtl_type="meQTL")

eqtl_expr <- read_tsv("results/current/biogrid_stringent/validation_expr/validation_all_eqtlgen.txt")
eqtl_tfa <- read_tsv("results/current/biogrid_stringent/validation_tfa/validation_all_eqtlgen.txt")
eqtl <- bind_rows(eqtl_expr,eqtl_tfa) %>%
  mutate(type=c(rep("expr", nrow(eqtl_expr)), rep("tfa", nrow(eqtl_tfa))),
         qtl_type="eQTL")

data <- bind_rows(meqtl, eqtl)

tfa_expr_plot <- data %>%
  ggplot(aes(x=reorder(graph_type, -cross_cohort_mcc, FUN=median), 
             y=cross_cohort_mcc, color=type)) + 
  geom_boxplot(position = "dodge") + 
  scb_binary + 
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  labs(title="",
       y="MCC",
       x="method",
       color = "measure:") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=0, vjust=0.5, angle=-45),
        legend.position = "bottom")
tfa_expr_plot

# ------------------------------------------------------------------------------
print("Figure 1 - Compile full plot")
# ------------------------------------------------------------------------------
ggarrange(simulation_mcc,
          tfa_expr_plot,
          ncol = 2, labels = c("A", "B"))

# ------------------------------------------------------------------------------
print("Done.\nSessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
