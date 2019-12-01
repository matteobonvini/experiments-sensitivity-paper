##################################################################
## Generate plot of the bounds for data of Connors et al (1996) ##
##################################################################
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("plot_theme.R")
library(RColorBrewer)

options(repr.plot.width = 10, repr.plot.height = 10)

dat <- read.csv("./results/data analysis/rhc_bounds.csv")

mm <- "x"
dat_mm <- dat[dat$model == mm & dat$delta == 1 & dat$epsilon <= 0.5, ]
xaxis_breaks <- c(0, dat_mm$eps_zero[1], 0.20, 0.30, 0.50)
unconfound_label <- as.character(round(dat_mm$lb[dat_mm$epsilon == 0][1] * 100, 
                                       1))
delta <- 1
ylims <- c(-0.45, 0.29)
eps0_label <- as.character(round(dat_mm$eps_zero[1] * 100, 1))
ytitle <- ifelse(mm == "x", "Difference in survival (%)", "")
xaxis_labs <- c(0, eps0_label, 20, 30, 50)

if(mm == "x") { colors <- brewer.pal(4, "YlOrRd") }
if(mm == "xa") { colors <- brewer.pal(4, "Blues") }

p <- ggplot() + 
  
  geom_ribbon(data = dat_mm, aes(x = epsilon, ymin = ci_lb_lo_unif, 
                                 ymax = ci_ub_hi_unif), 
              fill = colors[4], alpha = 1) +
  
  geom_ribbon(data = dat_mm, aes(x = epsilon, ymin = ci_im04_lo, 
                                 ymax = ci_im04_hi), fill = colors[2], 
              alpha = 1) +
  
  geom_ribbon(data = dat_mm, aes(x = epsilon, ymin = lb, ymax = ub), 
              fill = colors[1], alpha = 1) +
  
  geom_line(data = dat_mm, aes(x = epsilon, y = lb), color = "black", 
            size = 1.1) +
  
  geom_line(data = dat_mm, aes(x = epsilon, y = ub), color = "black",  
            size = 1.1) +
  
  geom_segment(data = dat_mm[1, ], size = 0.8, colour = "black", alpha = 0.5,
               aes(x = eps_zero, xend = eps_zero, y = ylims[1],  yend = 0)) +
  
  geom_segment(data = dat_mm[1, ], size = 0.8, colour = "black", alpha = 0.5,
               aes(x = 0, xend = eps_zero, y = 0, yend = 0)) +
  
  scale_x_continuous(breaks = xaxis_breaks,
                     expand = c(0, 0), limits = c(0, 0.502),
                     labels = xaxis_labs) +
  scale_y_continuous(breaks = c(-0.4, -0.2, dat_mm$lb[dat_mm$epsilon == 0][1],
                                0, 0.2),
                     expand = c(0,0),
                     limits = ylims,
                     labels = c(-40, -20, unconfound_label, 0, 20)) +
  
  scale_colour_manual(name = "model", values = colors) + 
  
  scale_fill_manual(name = "model", values = colors) + 
  
  labs(x = "% of Confounded Units", y = ytitle, colour = "", fill = "") +
  our_theme + 
  theme(legend.position = "none")

p

ggsave(filename = paste0("./results/data analysis/rhc_bounds_epsilon_", mm, 
                         ".pdf"), p)

## Superimposing the two plots
p2 <- ggplot() + 
  geom_ribbon(data = dat[dat$model == "xa", ], 
              aes(x = epsilon, ymin = ci_lo, ymax = ci_hi), fill = colors["xa"], 
              alpha = 0.25) +
  geom_ribbon(data = dat[dat$model == "xa", ], 
              aes(x=epsilon, ymin = ci_lo_im04, ymax = ci_hi_im04),
              fill = colors["xa"], alpha = 0.40) +
  geom_ribbon(data = dat[dat$model == "x", ], 
              aes(x = epsilon, ymin = lb, ymax = ub), 
              fill = colors["x"], alpha = 1) +
  geom_ribbon(aes(x = dat$epsilon[dat$model == "x"], 
                  ymin = dat$ub[dat$model == "x"], 
                  ymax = dat$ub[dat$model == "xa"]), fill = colors["xa"], 
              alpha = 1) +
  geom_ribbon(aes(x = dat$epsilon[dat$model == "x"], 
                  ymin = dat$lb[dat$model == "xa"], 
                  ymax = dat$lb[dat$model == "x"]), fill = colors["xa"], 
              alpha = 1) +
  geom_line(data = dat, aes(x = epsilon, y = lb, linetype = model), 
            color = "black", size = 1.1) +
  geom_line(data = dat, aes(x = epsilon, y = ub, linetype = model), 
            color = "black", size = 1.1) +
  geom_segment(data = dat[dat$model == "x", ][1, ], 
               aes(x = eps_zero, xend = eps_zero, y = -0.28, yend = 0), 
               size = 0.8, colour = "red", alpha = 0.75) +
  geom_segment(data = dat[dat$model == "xa", ][1, ], 
               aes(x = eps_zero, xend = eps_zero, y = -0.28, yend = 0), 
               size = 0.8, colour = "blue", alpha = 0.75) +
  geom_segment(data = dat[dat$model=="x", ][1, ], 
               aes(x = 0, xend = eps_zero, y = 0, yend = 0),
               size = 0.8, colour = "black", alpha = 0.5) +
  scale_x_continuous(breaks = c(0, dat[dat$model == "xa", ]$eps_zero[1], 
                              dat[dat$model == "x", ]$eps_zero[1], 0.2),
                     expand = c(0, 0), limits = c(0, 0.2),
                     labels = c("0", eps0_xa_label, eps0_x_label, "20")) +
  scale_y_continuous(breaks = c(-0.2, -0.1, dat$lb[dat$epsilon == 0][1], 
                                0, 0.1), expand = c(0,0), 
                     limits = c(-0.28, 0.20), 
                     labels = c("-20", "-10", unconfound_label, "0", "10")) +
  scale_colour_manual(name = "model", values = colors) + 
  scale_fill_manual(name = "model", values = colors) + 
  labs(x = "% of Confounded Units", y = "Difference in survival (%)", 
       colour = "", fill = "") +
  our_theme + 
  theme(legend.position = "none")

ggsave(filename = "./results/data analysis/rhc_bounds_epsilon.pdf", p2)

## Plot eps_delta
dat <- read.csv("./results/data analysis/rhc_bounds.csv")
mm <- "xa"
dat_mm <- dat[dat$model == mm & dat$epsilon <= 0.5, ]
dat_mm$delta <- as.factor(dat_mm$delta)
ytitle <- ifelse(mm == "x", "Difference in survival (%)", "")
deltas <- sort(unique(dat_mm$delta))
names(colors) <- deltas

eps0_label1 <- as.character(round(dat_mm[dat_mm$delta == deltas[1], ]$eps_zero[1] * 100))
eps0_label2 <- as.character(round(dat_mm[dat_mm$delta == deltas[2], ]$eps_zero[1] * 100))
eps0_label3 <- as.character(round(dat_mm[dat_mm$delta == deltas[3], ]$eps_zero[1] * 100))
eps0_label4 <- as.character(round(dat_mm[dat_mm$delta == deltas[4], ]$eps_zero[1] * 100))

if(mm == "x") {
  colors <- rev(brewer.pal(4, "YlOrRd"))
  newcol <- "blue"
  xaxis_col <- c("black", newcol, newcol, newcol, "black", newcol, "black")
  xaxis_bold <- ifelse(xaxis_col == "black", "plain" , "bold")
  ylims <- c(-0.45, 0.29)
  xaxis_labs <- c(0, eps0_label4, eps0_label3, eps0_label2, 25,
                  eps0_label1, 50)
  xaxis_breaks <- c(0, 
                    dat_mm[dat_mm$delta == deltas[4], ]$eps_zero[1],
                    dat_mm[dat_mm$delta == deltas[3], ]$eps_zero[1],
                    dat_mm[dat_mm$delta == deltas[2], ]$eps_zero[1],
                    0.25,
                    dat_mm[dat_mm$delta == deltas[1], ]$eps_zero[1],
                    0.5)
} 
if(mm == "xa")  {
  colors <- rev(brewer.pal(4, "Blues"))
  newcol <- "red"
  xaxis_col <- c("black", newcol, newcol, newcol, "black", newcol, "black")
  xaxis_bold <- ifelse(xaxis_col == "black", "plain" , "bold")
  ylims <- c(-0.45, 0.29)
  xaxis_labs <- c(0, eps0_label4, eps0_label3, eps0_label2, 25,
                  eps0_label1, 50)
  xaxis_breaks <- c(0, 
                    dat_mm[dat_mm$delta == deltas[4], ]$eps_zero[1],
                    dat_mm[dat_mm$delta == deltas[3], ]$eps_zero[1],
                    dat_mm[dat_mm$delta == deltas[2], ]$eps_zero[1],
                    0.25,
                    dat_mm[dat_mm$delta == deltas[1], ]$eps_zero[1],
                    0.5)
}

p3 <- ggplot() + 
  
  geom_ribbon(data = dat_mm[dat_mm$delta == deltas[1], ], 
              aes(x = epsilon, ymin = lb, ymax = ub), 
              fill = colors[deltas[1]], alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[2]], 
                  ymin = dat_mm$ub[dat_mm$delta == deltas[1]], 
                  ymax = dat_mm$ub[dat_mm$delta == deltas[2]]), 
              fill = colors[deltas[2]], 
              alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[2]], 
                  ymin = dat_mm$lb[dat_mm$delta == deltas[2]], 
                  ymax = dat_mm$lb[dat_mm$delta == deltas[1]]), 
              fill = colors[deltas[2]], 
              alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[3]], 
                  ymin = dat_mm$ub[dat_mm$delta == deltas[2]], 
                  ymax = dat_mm$ub[dat_mm$delta == deltas[3]]), 
              fill = colors[deltas[3]], 
              alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[3]], 
                  ymin = dat_mm$lb[dat_mm$delta == deltas[3]], 
                  ymax = dat_mm$lb[dat_mm$delta == deltas[2]]), 
              fill = colors[deltas[3]], 
              alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[4]], 
                  ymin = dat_mm$ub[dat_mm$delta == deltas[3]], 
                  ymax = dat_mm$ub[dat_mm$delta == deltas[4]]), 
              fill = colors[deltas[4]], 
              alpha = 1) +
  
  geom_ribbon(aes(x = dat_mm$epsilon[dat_mm$delta == deltas[4]], 
                  ymin = dat_mm$lb[dat_mm$delta == deltas[4]], 
                  ymax = dat_mm$lb[dat_mm$delta == deltas[3]]), 
              fill = colors[deltas[4]], 
              alpha = 1) +
  
  geom_line(data = dat_mm, aes(x = epsilon, y = lb, linetype = delta), 
            size = 0.8, color = "black") +
  
  geom_line(data = dat_mm, aes(x = epsilon, y = ub, linetype = delta), 
            size = 0.8, color = "black") +

  geom_segment(data = dat_mm[dat_mm$delta == deltas[1], ], size = 0.8,
               colour = newcol, alpha = 0.5,
               aes(x = 0, xend = 0.50, y = 0, yend = 0)) +

  scale_x_continuous(breaks = xaxis_breaks,
                     expand = c(0, 0), limits = c(0, 0.502),
                     labels = xaxis_labs) +
  
  scale_y_continuous(breaks = c(-0.4, -0.2, dat_mm$lb[dat_mm$epsilon == 0][1],
                                0, 0.2),
                     expand = c(0,0),
                     limits = ylims,
                     labels = c(-40, -20, unconfound_label, 0, 20)) +
  
  labs(x = "% of Confounded Units", y = ytitle, colour = "delta", fill = "") +
  our_theme +
  theme(axis.text.x = element_text(colour = xaxis_col, face = xaxis_bold)) +
  theme(legend.position = "none")
p3
ggsave(filename = paste0("./results/data analysis/rhc_bounds_epsilon_delta_", mm, 
                         ".pdf"), p3)
  
