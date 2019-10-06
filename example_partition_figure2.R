rm(list = ls())

set.seed(1000)

library(tidyverse)
library(RColorBrewer)
library(ggalt)
library(sensAteBounds)
library(truncnorm)

plotsize <- function(x,y) options(repr.plot.width=x, repr.plot.height=y)
plotsize(10, 10)
colors <- brewer.pal(9, "Set1")
colors <- colors[c(1,2,7,9)]

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plot_theme.R")
source("true_regression_functions_simulation.R")

truth <- readRDS("./data/truth_simulation.RData")
df <- gen_data(n=5000)
gx_vals <- gx(df$x1, df$x2)

# Select epsilon points for example in Figure 2
eps_seq <- c(0.01, 0.05, 0.15)
bounds <- round(truth[truth[, "eps"] %in% eps_seq, ], 2)

# Generate plots in Figure 2 for each epsilon
for(eps in eps_seq) {
  
  shatl <- I(gx_vals  <= quantile(gx_vals , eps))
  shatu <- I(gx_vals > quantile(gx_vals, 1-eps))
  df$group <- factor("none", levels = c("none", "lower bound", "upper bound"))
  df$group[shatl] <- "lower bound"
  df$group[shatu] <- "upper bound"
  
  ggplot(data = df, aes(x1, x2)) + 
    geom_point(aes(colour = group)) + 
    scale_colour_manual(values=c("grey", colors[1], colors[2]), name="colour") + 
    labs(x = "Covariate 1", y = "Covariate 2", colour = "") +
    scale_x_continuous(breaks=c(-2, 0, 2),
                       labels=c("-2", "0", "2")) +
    scale_y_continuous(breaks=c(-2, 0, 2),
                       labels=c("-2", "0", "2")) +
    our_theme + 
    theme(legend.position = "none")
  ggsave(filename = paste0("./results/simulation/clusters_", 100 * eps, ".pdf"))
}
