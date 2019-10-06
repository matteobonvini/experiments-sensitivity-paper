rm(list = ls())
library(RColorBrewer)

setwd("C:/Users/matte/Dropbox/Causal questions shared with Matteo/Sensitivity Analysis/experiments")
source("./plot_theme.R")

plotsize <- function(x,y) options(repr.plot.width=x, repr.plot.height=y)
plotsize(10, 10)
colors <- brewer.pal(8, "Set2")
colors <- c("#E69F00", "#56B4E9")

dat <- read.csv("./results/data analysis/rhc_bounds_temp.csv")
names(colors) <- c("x", "xa")

ggplot() + 
  geom_ribbon(data = dat[dat$model == "xa", ], aes(x=epsilon, ymin=ci_lo, ymax=ci_hi), 
              fill = colors["xa"], alpha=0.25) +
  geom_ribbon(data =dat[dat$model == "xa", ], aes(x=epsilon, ymin=ci_lo_ptw_im04, ymax=ci_hi_ptw_im04),
              fill = colors["xa"], alpha=0.40) +
  geom_ribbon(data = dat[dat$model == "x", ], aes(x=epsilon, ymin=lb, ymax=ub), 
              fill = colors["x"], alpha = 1) +
  geom_ribbon(aes(x=dat$epsilon[dat$model == "x"], ymin=dat$ub[dat$model == "x"], ymax=dat$ub[dat$model == "xa"]), 
              fill = colors["xa"], alpha = 1) +
  geom_ribbon(aes(x=dat$epsilon[dat$model == "x"], ymin=dat$lb[dat$model == "xa"], ymax=dat$lb[dat$model == "x"]), 
              fill = colors["xa"], alpha = 1) +
  geom_line(data = dat, aes(x=epsilon, y=lb, linetype=model), 
            color = "black", size=1.1) +
  geom_line(data = dat, aes(x=epsilon, y=ub, linetype=model), 
            color = "black", size=1.1) +
  geom_segment(data = dat[dat$model == "x", ][1, ], 
               aes(x=eps_zero, xend=eps_zero, y=-0.28, yend=0), 
               size=0.8, colour="red", alpha=0.75) +
  geom_segment(data = dat[dat$model == "xa", ][1, ], 
               aes(x=eps_zero, xend=eps_zero, y=-0.28, yend=0), 
               size=0.8, colour="blue", alpha=0.75) +
  geom_segment(data = dat[dat$model=="x", ][1, ], 
               aes(x=0, xend=eps_zero, y=0, yend=0),
               size=0.8, colour="black", alpha=0.5) +
  scale_x_continuous(breaks=c(0, dat[dat$model=="xa", ]$eps_zero[1], 
                              dat[dat$model=="x", ]$eps_zero[1], 0.2),
                     expand = c(0,0), limits = c(0, 0.2),
                     labels=c("0", as.character(round(dat[dat$model=="xa", ]$eps_zero[1]*100)), 
                              as.character(round(dat[dat$model=="x", ]$eps_zero[1]*100)), "20")) +
  scale_y_continuous(breaks=c(-0.2, -0.1, dat$lb[dat$epsilon==0][1], 0, 0.1),
                     expand = c(0,0), limits = c(-0.28, 0.20),
                     labels=c("-20", "-10", as.character(round(dat$lb[dat$epsilon==0][1]*100, 1)), "0", "10")) +
  scale_colour_manual(name = "model", values=colors) + 
  scale_fill_manual(name = "model", values=colors) + 
  labs(x="% of Confounded Units", y="Difference in survival (%)", 
       colour="", fill="") +
  our_theme + 
  theme(legend.position = "none")

ggsave(filename = "./results/data analysis/rhc_bounds_epsilon.pdf")
