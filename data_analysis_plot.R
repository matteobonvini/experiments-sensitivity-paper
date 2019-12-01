##################################################################
## Generate plot of the bounds for data of Connors et al (1996) ##
##################################################################
rm(list = ls())

library(RColorBrewer)


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("plot_theme.R")

dat <- read.csv("./results/data analysis/rhc_bounds.csv")
datx <- dat[dat$model == "x" & dat$delta == 1 & dat$epsilon <= 0.5, ]
datxa <- dat[dat$model == "xa" & dat$delta == 1 & dat$epsilon <= 0.5, ]

ylims <- c(-0.45, 0.32)
models <- c("x", "xa")

colors_x <- brewer.pal(4, "YlOrRd")
colors_xa <- brewer.pal(4, "Blues")

eps0_label_fn <- function(dat, delta) {
  # Return value for estimate of eps0 for a given delta, dat must be data.frame
  # that results from get_bound() function. 
  out <- round(dat[dat$delta == delta, ]$eps_zero[1], 2)
  return(out)
}

unconfound_val <- round(dat$lb[dat$epsilon == 0][1], 3)
unconfound_label <- as.character(unconfound_val * 100)

eps0_xa_label <- as.character(eps0_label_fn(datxa, 1) * 100)
eps0_x_label <- as.character(eps0_label_fn(datx, 1) * 100)

for(mm in models) {
  ## Generate plots for Figure 1, eps in [0, 0.5], delta = 1.
  dat_mm <- dat[dat$model == mm & dat$delta == 1 & dat$epsilon <= 0.5, ]
  xaxis_breaks <- c(0, dat_mm$eps_zero[1], 0.20, 0.30, 0.50)

  ytitle <- ifelse(mm == "x", "Difference in survival (%)", "")

  if(mm == "x") { 
    colors <- colors_x 
    eps0_label <- eps0_x_label
    }
  if(mm == "xa") { 
    colors <- colors_xa
    eps0_label <- eps0_xa_label
  }
  
  xaxis_labs <- c(0, eps0_label, 20, 30, 50)
  
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
    scale_y_continuous(breaks = c(-0.4, -0.2, unconfound_val, 0, 0.2),
                       expand = c(0,0), limits = ylims,
                       labels = c(-40, -20, unconfound_label, 0, 20)) +
    
    scale_colour_manual(name = "model", values = colors) + 
    
    scale_fill_manual(name = "model", values = colors) + 
    
    labs(x = "% of Confounded Units", y = ytitle, colour = "", fill = "") +
    our_theme + 
    theme(legend.position = "none")
  
  ggsave(filename = paste0("./results/data analysis/rhc_bounds_epsilon_", mm, 
                           ".pdf"), p, height = 3, width = 3)
}

## Superimposing the two plots
p2 <- ggplot() + 
  
  geom_ribbon(data = datxa, 
              aes(x = epsilon, ymin = lb, ymax = ub), fill = colors_xa[2], 
              alpha = 1) +
  
  geom_ribbon(data = datx, 
              aes(x = epsilon, ymin = lb, ymax = ub), fill = colors_x[2], 
              alpha = 1) +
  
  geom_line(data = datx, aes(x = epsilon, y = lb), color = "black", 
            size = 1.1) +
  
  geom_line(data = datxa, aes(x = epsilon, y = lb), color = "black", 
            linetype = "dashed", size = 1.1) +
  
  geom_line(data = datx, aes(x = epsilon, y = ub), color = "black", 
            size = 1.1) +
  
  geom_line(data = datxa, aes(x = epsilon, y = ub), color = "black", 
            linetype = "dashed", size = 1.1) +
  
  geom_segment(data = datx[1, ], 
               aes(x = eps_zero, xend = eps_zero, y = ylims[1], yend = 0), 
               size = 0.8, colour = "red", alpha = 0.75) +
  
  geom_segment(data = datxa[1, ], 
               aes(x = eps_zero, xend = eps_zero, y = ylims[1], yend = 0), 
               size = 0.8, colour = "blue", alpha = 0.75) +
  
  geom_segment(data = datx[1, ], 
               aes(x = 0, xend = eps_zero, y = 0, yend = 0),
               size = 0.8, colour = "black", alpha = 0.5) +
  
  scale_x_continuous(breaks = c(0, datxa$eps_zero[1],
                                datx$eps_zero[1], 0.2,0.5),
                     expand = c(0, 0), limits = c(0, 0.5),
                     labels = c(0, eps0_xa_label, eps0_x_label, 20, 50)) +
  
  scale_y_continuous(breaks = c(-0.2, -0.1, unconfound_val, 0, 0.1), 
                     expand = c(0,0), limits = ylims, 
                     labels = c("-20", "-10", unconfound_label, "0", "10")) +
  
  scale_colour_manual(name = "model", values = colors) + 
  
  scale_fill_manual(name = "model", values = colors) + 
  
  labs(x = "% of Confounded Units", y = "Difference in survival (%)", 
       colour = "", fill = "") +
  
  our_theme + 
  
  theme(legend.position = "none")

ggsave(filename = "./results/data analysis/rhc_bounds_epsilon.pdf", p2,
       height = 3, width = 3)

## Plot eps_delta
dat <- read.csv("./results/data analysis/rhc_bounds.csv")
models <- c("x", "xa")
ylims <- c(-0.45, 0.32)

for(mm in models) {
  ## Generate plots for appendix figure, eps in [0, 0.5] with different deltas.
  dat_mm <- dat[dat$model == mm & dat$epsilon <= 0.5, ]
  dat_mm$delta <- as.factor(dat_mm$delta)
  
  # Since figure for "xa" is plotted to the right of that for "x", we avoid 
  # repeating ylab.
  ytitle <- ifelse(mm == "x", "Difference in survival (%)", "")
  
  deltas <- sort(unique(dat_mm$delta))
  names(colors) <- deltas
  
  eps0_vals <- sapply(deltas, eps0_label_fn, dat = dat_mm)
  xaxis_breaks <- sort(c(0, eps0_vals[4], eps0_vals[3], eps0_vals[2], 0.25,
                       eps0_vals[1], 0.5))
  xaxis_labs <- as.character(xaxis_breaks * 100)
  
  if(mm == "x") {
    colors <- rev(brewer.pal(4, "YlOrRd"))
    newcol <- "blue"
    xaxis_col <- c("black", newcol, newcol, newcol, "black", newcol, "black")
    xaxis_bold <- ifelse(xaxis_col == "black", "plain" , "bold")
  } 
  if(mm == "xa")  {
    colors <- rev(brewer.pal(4, "Blues"))
    newcol <- "red"
    xaxis_col <- c("black", newcol, newcol, newcol, "black", newcol, "black")
    xaxis_bold <- ifelse(xaxis_col == "black", "plain" , "bold")
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
  filename <- paste0("./results/data analysis/rhc_bounds_epsilon_delta_", mm, 
                     ".pdf")
  ggsave(filename = filename, p3, height = 3, width = 3)
}

