# Design sensitivity
rm(list = ls())

set.seed(1000)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("simulation_true_regression_functions.R")

sim_fun <- function(alpha) {
  
  expect_y0 <- function(x1, x2, u){
    return(0.5)
  }
  
  mu0x <- function(x1, x2) {
    u0 <- expect_y0(x1, x2, 0) * (1-pi_xu(x1, 0)) / (1-pix(x1))
    u1 <- expect_y0(x1, x2, 1) * (1-pi_xu(x1, 1)) / (1-pix(x1))
    return((1-pu_x(x1)) * u0 + pu_x(x1) * u1)
  }
  
  mu0 <- adaptIntegrate(f=function(x) { 
    mu0x(x[1], x[2]) * densx(x[1], x[2]) 
  }, lowerLimit=c(xlb, xlb), upperLimit=c(xub, xub), absError=1e-15)$integral
  
  expect_y1 <- function(x1, x2, u){
    return( qbeta(ptruncnorm(x2, xlb, xub), alpha, 1))
  }
  
  mu1x <- function(x1, x2) {
    u0 <- expect_y1(x1, x2, 0) * pi_xu(x1, 0) / pix(x1)
    u1 <- expect_y1(x1, x2, 1) * pi_xu(x1, 1) / pix(x1)
    return((1-pu_x(x1)) * u0 + pu_x(x1) * u1)
  }
  
  mu1 <- adaptIntegrate(f=function(x) { 
    mu1x(x[1], x[2]) * densx(x[1], x[2]) 
  }, lowerLimit=c(xlb, xlb), upperLimit=c(xub, xub), absError=1e-15)$integral
  
  x1vals <- rtruncnorm(1e4, xlb, xub)
  x2vals <- rtruncnorm(1e4, xlb, xub)
  avals <- rbinom(1e4, 1, pix(x1vals))
  
  glb_x <- gx(x1vals, x2vals) - 1
  glb_xa <- gxa(x1vals, x2vals, avals) - 1
  
  quants_lb_x <- quantile(glb_x, eps_seq)
  quants_lb_xa <- quantile(glb_xa, eps_seq)
  
  g_term_lb_x <- get_g_term(g = glb_x, quants = quants_lb_x, upper = FALSE)
  g_term_lb_xa <- get_g_term(g = glb_xa, quants = quants_lb_xa, upper = FALSE)
  
  lb_x <- mu1 - mu0 + g_term_lb_x
  lb_xa <- mu1 - mu0 + g_term_lb_xa
  
  eps0_x <- eps_seq[which.min(abs(lb_x))]
  eps0_xa <- eps_seq[which.min(abs(lb_xa))]
  
  out <- array(c(mu1 - mu0, eps0_x, eps0_xa), dim = c(1, 3))
  
  return(out)
  
}

eps_seq <- seq(0, 1, 0.01)
nsim <- 100
alphas <- seq(1, 50, length.out = 100)

res <- matrix(NA, ncol = 3, nrow = length(alphas), 
              dimnames = list(NULL, c("tau", "eps_x", "eps_xa")))
res_sd <- matrix(NA, ncol = 3, nrow = length(alphas), 
              dimnames = list(NULL, c("tau", "eps_x", "eps_xa")))
for(i in 1:length(alphas)) {
  sim_res <- replicate(nsim, sim_fun(alphas[i]))
  res[i, ] <- apply(sim_res, c(1, 2), mean)
  res_sd[i, ] <- apply(sim_res, c(1, 2), sd) / sqrt(nsim)
  
  print(paste0("Progress = ", round(i / length(alphas) * 100, 2), " %"))
}

plot(x = alphas, y = res[, 1], type = "l", ylim = c(0, 1))
lines(x = alphas, y = res[, 2], col = "red")
lines(x = alphas, y = res[, 3], col = "blue")

saveRDS(res, file = "./results/simulation/design_sensitivity_res.RData")
library(ggplot2)
source("plot_theme.R")

label_curves <- c(bquote(tau), 
                  bquote(paste(tilde(epsilon), " under ",
                                   X, "-model")),
                  bquote(paste(tilde(epsilon), " under ",
                                   XA, "-model")))

p <- ggplot(NULL, aes(x = alphas)) +
  geom_line(aes(y = res[, "tau"], colour = "tau", linetype = "tau"), 
            size = 2.5) +
  geom_line(aes(y = res[, "eps_x"], colour = "x", linetype = "x"), 
            size = 2.5) +
  geom_line(aes(y = res[, "eps_xa"], colour = "xa", linetype = "xa"), 
            size = 2.5) + 
  scale_x_continuous(expand = c(0.05, 0.05), 
                     breaks = c(1, 25, 50)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = expression(alpha), y = "", color = "", fill = "") +
  scale_colour_manual(values = c("tau" = "black", "x" = "red",
                                 "xa" = "blue"), name = "", 
                      labels = label_curves) +
  scale_linetype_manual(values=c("tau" = "solid",
                                 "x" = "dashed", 
                                 "xa" = "dotted"),
                        name = "",
                        labels = label_curves) +
  our_theme + 
  theme(legend.position = c(0.50, 0.2), 
        legend.text = element_text(colour = "black", size=50, face="bold"),
        legend.direction = "vertical",
        axis.title.x = element_text(size = 50),
        axis.text.x = element_text(size = 50),
        axis.text.y = element_text(size = 50))
p
ggsave(filename = paste0("./results/simulation/design_sensitivity.pdf"), p, 
       height = 10, width = 10)
