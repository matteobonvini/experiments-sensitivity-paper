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
source("simulation_true_regression_functions.R")

truth <- readRDS("./data/truth_simulation.RData")
df <- gen_data(n=5000)
gx_vals <- gx(df$x1, df$x2)

# Select epsilon points for example in Figure 2
eps_seq <- c(0.01, 0.05, 0.15)
bounds <- round(truth[truth[, "eps"] %in% eps_seq, ], 2)

# Color units the yield lower and upper bounds from a model as in the
# simulation section.
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

## Generate example of a partition
rm(list = ls())

plotsize <- function(x,y) options(repr.plot.width=x, repr.plot.height=y)
plotsize(10, 10)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("plot_theme.R")

n <- 100
x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1)

plot1 <- ggplot() + 
  geom_point(aes(x = x1, y = x2), size = 2.5) +
  labs(x = "Covariate 1", y = "Covariate 2", colour = "") +
  scale_x_continuous(breaks=c(0, 0.5, 1)) +
  scale_y_continuous(breaks=c(0, 0.5, 1)) +
  our_theme
ggsave(plot1, filename = paste0("../draft/figures/ex_partition.jpeg"))

# Then I manually draw the partitions in ppt, since this is just for an example. 
