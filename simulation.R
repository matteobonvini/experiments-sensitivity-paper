rm(list = ls())

set.seed(1000)

library(pbapply)
library(devtools)
devtools::install("C:/Users/matte/Desktop/sensAteBounds")
setwd("C:/Users/matte/Desktop/sensAteBounds")
devtools::test()
library(sensAteBounds)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("true_regression_functions_simulation.R")

## Load true values ##
truth <- readRDS("./data/truth_simulation.RData")
eps_seq <- attributes(truth)$eps_seq
eps0_seq <- attributes(truth)$eps0_seq

# SL library
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam")
nsim <- 500
n <- c(10000, 5000, 1000, 500)
alpha <- 0.05

sim_fn <- function(n) {
  # Function to simulate data and estimate bounds under S \ind (Y, A) | X  
  df <- gen_data(n)
  y <- df$y
  a <- df$a 
  x <- df[, c("x1", "x2")]
  
  # estimate the nuisance regression functios
  nuis_fns <- do_crossfit(y = y, a = a, x = x, nsplits = 2, outfam = binomial(),
                          treatfam = binomial(), sl.lib = sl.lib)
  
  res <- get_bound(y = y, a = a, x = x, outfam = NULL, treatfam = NULL, 
                   model = "x", eps = eps0_seq, delta = 1, nsplits = NULL, 
                   do_mult_boot = FALSE, do_eps_zero = TRUE, 
                   nuis_fns = nuis_fns, alpha = alpha)
  bounds <- res$bounds
  eps_zero <- res$eps_zero
  # select subset of esp0_seq where to evaluate the bounds curves
  idx <- which(round(bounds$eps, 10) %in% round(eps_seq, 10))
  
  lb <- bounds$lb[idx]
  ub <- bounds$ub[idx]
  
  idx1 <- which(round(res$phibar_lb$eps, 10) %in% round(eps_seq, 10)) 
  phibar_lb <- res$phibar_lb[idx1, ]
  phibar_ub <- res$phibar_ub[idx1, ]
  
  idx2 <- which(round(res$lb_var$eps, 10) %in% round(eps_seq, 10)) 
  lb_var <- res$lb_var[idx2, ]
  ub_var <- res$ub_var[idx2, ]
  
  lb_estbar <- summarize_by(phibar_lb, "mean", delta, eps)
  ub_estbar <- summarize_by(phibar_ub, "mean", delta, eps)
  
  calpha_lb <- get_multboot(n=n, psihat=lb_estbar$ifvals,
                            sigmahat=lb_var$ifvals,
                            ifvals=phibar_lb$ifvals, alpha=alpha/2, B=1e4)
  calpha_ub <- get_multboot(n=n, psihat=-ub_estbar$ifvals,
                            sigmahat=ub_var$ifvals,
                            ifvals=-phibar_ub$ifvals, alpha=alpha/2, B=1e4)
  
  ci_lb <- get_ci(lb, sqrt(lb_var$ifvals/n), calpha_lb)
  ci_ub <- get_ci(ub, sqrt(ub_var$ifvals/n), calpha_ub)
  
  cnames <- c("eps", "n", "lb", "ub", "eps_zero", "ci_lo", "ci_lb_hi",
              "ci_ub_lo", "ci_hi", "eps_zero_lo", "eps_zero_hi")
  nidx <- length(bounds$eps[idx])
  out <- c(bounds$eps[idx], rep(n, nidx), lb, ub, rep(eps_zero$eps_zero, nidx),
           ci_lb[, 1], ci_lb[, 2], ci_ub[, 1], ci_ub[, 2], 
           rep(eps_zero$eps_zero_lo, nidx), rep(eps_zero$eps_zero_hi, nidx))
  out <- matrix(out, ncol = length(cnames), nrow = nidx, 
                dimnames = list(NULL, cnames))
  return(out)
}

# Store the results for Table 1 in Section 5
res_cnames <- c("n", "bias_lb", "bias_ub", "bias_eps0", "rmse_lb", "rmse_ub", 
                "cvg_reg", "cvg_eps0")
res <- matrix(NA, ncol = length(res_cnames), nrow = length(n), 
              dimnames = list(NULL, res_cnames))

for (i in 1:length(n)) {
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn(n[i]))
  res[i, "n"] <- n[i]
  
  res[i, "bias_lb"] <- bias(sims[, "lb", ], truth[, "lb"])
  res[i, "bias_ub"] <- bias(sims[, "ub", ], truth[, "ub"])
  res[i, "bias_eps0"] <- bias(sims[, "eps_zero", ], truth[, "eps_zero"])
  
  res[i, "rmse_lb"] <- rmse(sims[, "lb", ], truth[, "lb"], n[i])
  res[i, "rmse_ub"] <- rmse(sims[, "ub", ], truth[, "ub"], n[i])
  
  res[i, "cvg_reg"] <- coverage(sims[, "ci_lo", ], sims[, "ci_hi", ], 
                                truth[, "lb"], truth[, "ub"])
  res[i, "cvg_eps0"] <- coverage(sims[, "eps_zero_lo", ], 
                                 sims[, "eps_zero_hi", ], 
                                 truth[, "eps_zero"], truth[, "eps_zero"])
  
  print(res[1:i, ])
  
  saveRDS(sims, file=paste0("./results/simulation/sims", n[i], ".RData"))
  saveRDS(res, file=paste0("./results/simulation/sim_res", n[i], ".RData"))
}

saveRDS(res, file=paste0("./results/simulation/sim_res.RData"))
