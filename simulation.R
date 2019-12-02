####################################
## Run simulation as in Section 5 ##
####################################
rm(list = ls())

library(sensitivitypuc)

set.seed(1000)

# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)
# devtools::install_github("matteobonvini/sensitivitypuc")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("simulation_true_regression_functions.R")

## Load true values ##
truth <- readRDS("./data/truth_simulation.RData")
# Sequence to evaluate the curves
eps_seq <- attributes(truth)$eps_seq
# Sequence (longer) to evaluate epsilon_zero
eps0_seq <- attributes(truth)$eps0_seq

# SuperLearner library
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam")
nsim <- 10
n <- rev(c(500, 1000, 5000, 10000))
n <- 500
alpha <- 0.05

sim_fn <- function(n) {
  # Function to simulate data and estimate bounds under S \ind (Y, A) | X  
  df <- gen_data(n)
  y <- df$y
  a <- df$a 
  x <- df[, c("x1", "x2")]
  
  res <- get_bound(y = y, a = a, x = x, outfam = binomial(), model = "x",
                   ymin = 0, ymax = 1, treatfam = binomial(), eps = eps0_seq, 
                   delta = 1, nsplits = 5, do_mult_boot = FALSE, 
                   do_eps_zero = TRUE, alpha = alpha, sl.lib = sl.lib, 
                   do_parallel = TRUE, ncluster = 3, do_rearrange = FALSE,
                   show_progress = FALSE)
  
  bounds <- res$bounds[, , 1]
  eps_zero <- res$eps_zero
  
  # select subset of esp0_seq where to evaluate the bounds curves
  idx <- which(round(eps0_seq, 10) %in% round(eps_seq, 10))
  
  lb <- bounds[idx, "lb"]
  ub <- bounds[idx, "ub"]
  
  phibar_l <- res$phibar_lb[, idx, ]
  phibar_u <- res$phibar_ub[, idx, ]
  
  var_l <- res$var_lb[idx, 1]
  var_u <- res$var_ub[idx, 1]
  
  lb_estbar <- apply(phibar_l, 2, mean)
  ub_estbar <- apply(phibar_u, 2, mean)
  
  calpha_lb <- do_multboot(n = n, psihat = lb_estbar, sigmahat = var_l, 
                            ifvals = phibar_l, alpha = alpha/2, B = 10000)
  calpha_ub <- do_multboot(n = n, psihat = -ub_estbar, sigmahat = var_u,
                            ifvals = -phibar_u, alpha = alpha/2, B = 10000)
  
  
  ci_lb <- get_ci(lb, sqrt(var_l/n), calpha_lb)
  ci_ub <- get_ci(ub, sqrt(var_u/n), calpha_ub)
  
  # apply rearrangement of Chernozhukov (2009)
  lb <- sort(lb, decreasing = TRUE)
  ub <- sort(ub)
  ci_lb <- apply(ci_lb, 2, sort, decreasing = TRUE)
  ci_ub <- apply(ci_ub, 2, sort)
  
  cnames <- c("eps", "n", "lb", "ub", "eps_zero", "ci_lo", "ci_lb_hi",
              "ci_ub_lo", "ci_hi", "eps_zero_lo", "eps_zero_hi")
  n_row <- length(eps_seq)
  n_col <- length(cnames) 
  out <- c(eps_seq, rep(n, n_row), lb, ub, rep(eps_zero$est, n_row),
           ci_lb[, 1], ci_lb[, 2], ci_ub[, 1], ci_ub[, 2], 
           rep(eps_zero$ci_lo, n_row), rep(eps_zero$ci_hi, n_row))
  out <- matrix(out, ncol = n_col, nrow = n_row, dimnames = list(NULL, cnames))
  
  return(out)
}

# Store the results for Table 1 in paper
res_cnames <- c("n", "bias_lb", "bias_ub", "bias_eps0", "rmse_lb", "rmse_ub", 
                "rmse_eps0", "cvg_reg", "cvg_eps0")
res <- matrix(NA, ncol = length(res_cnames), nrow = length(n), 
              dimnames = list(NULL, res_cnames))

for (i in 1:length(n)) {
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn(n[i]))
  res[i, "n"] <- n[i]
  
  res[i, "bias_lb"] <- bias(sims[, "lb", ], truth[, "lb"])
  res[i, "bias_ub"] <- bias(sims[, "ub", ], truth[, "ub"])
  res[i, "bias_eps0"] <- bias(sims[1, "eps_zero", ], truth[1, "eps_zero"])
  
  res[i, "rmse_lb"] <- rmse(sims[, "lb", ], truth[, "lb"], n[i])
  res[i, "rmse_ub"] <- rmse(sims[, "ub", ], truth[, "ub"], n[i])
  res[i, "rmse_eps0"] <- rmse(sims[1, "eps_zero", ], truth[1, "eps_zero"], n[i])
  
  res[i, "cvg_reg"] <- coverage(sims[, "ci_lo", ], sims[, "ci_hi", ], 
                                truth[, "lb"], truth[, "ub"])
  res[i, "cvg_eps0"] <- coverage(sims[1, "eps_zero_lo", ], 
                                 sims[1, "eps_zero_hi", ], 
                                 truth[1, "eps_zero"], truth[1, "eps_zero"])

  print(res[1:i, ])
  
  # saveRDS(sims, file=paste0("./results/simulation/sims", n[i], ".RData"))
  # saveRDS(res, file=paste0("./results/simulation/sim_res", n[i], ".RData"))
}

# saveRDS(res, file=paste0("./results/simulation/sim_res.RData"))

# print(round(readRDS(paste0("./results/simulation/sim_res.RData")), 2))
