####################################
## Run simulation as in Section 5 ##
####################################
rm(list = ls())
library(devtools) 
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # if warnings are causing erros and they are not important one may use this
devtools::install_github("matteobonvini/sensitivitypuc") 
library(sensitivitypuc)

set.seed(1000)

library(sensitivitypuc)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("simulation_true_regression_functions.R")

## Load true values ##
truth <- readRDS("./data/truth_simulation.RData")
# Sequence to evaluate the curves
eps_seq <- attributes(truth)$eps_seq
# Sequence (longer) to evaluate epsilon_zero
eps0_seq <- attributes(truth)$eps0_seq

# SuperLearner library
sl.lib <- c("SL.mean", "SL.speedlm", "SL.speedglm", "SL.glm.interaction", 
            "SL.gam")
nsim <- 500
n <- rev(c(500, 1000, 5000, 10000))
nsplits <- c(5, 5, 5, 5)
alpha <- 0.05

sim_fn <- function(n, nsplits) {
  # Function to simulate data and estimate bounds under S \ind (Y, A) | X  
  df <- gen_data(n)
  y <- df$y
  a <- df$a 
  x <- df[, c("x1", "x2")]
  ymin <- 0
  ymax <- 1
  delta <- 1
  outfam <- treatfam <- binomial()
  model <- "x"
  show_progress <- FALSE
  do_rearrange <- TRUE
  res_zero <- get_bound(y = y, a = a, x = x, model = model, ymin = ymin, 
                        ymax = ymax, outfam = outfam, treatfam = treatfam, 
                        eps = eps0_seq, delta = delta, nsplits = nsplits, 
                        do_mult_boot = FALSE, do_eps_zero = TRUE, alpha = alpha, 
                        sl.lib = sl.lib, do_parallel = TRUE, ncluster = 3, 
                        do_rearrange = do_rearrange, 
                        show_progress = show_progress)
  
  nuis_fns <- res_zero$nuis_fn
  bounds <- res_zero$bounds[, , 1]
  eps_zero <- res_zero$eps_zero
  
  res <- get_bound(y = y, a = a, x = x, model = model, ymin = ymin, 
                   ymax = ymax, outfam = NULL, treatfam = NULL, 
                   eps = eps_seq, delta = delta, nsplits = NULL, 
                   do_mult_boot = TRUE, B = 10000, do_eps_zero = FALSE, 
                   alpha = alpha, sl.lib = NULL, ncluster = NULL,
                   do_rearrange = do_rearrange, do_parallel = FALSE,
                   show_progress = show_progress, nuis_fns = nuis_fns)

  lb <- res$bounds[, "lb", ]
  ub <- res$bounds[, "ub", ]
  
  ci_lb_lo <- res$bounds[, "ci_lb_lo_unif", ]
  ci_ub_hi <- res$bounds[, "ci_ub_hi_unif", ]
  
  cnames <- c("eps", "n", "lb", "ub", "eps_zero", "ci_lb_lo", "ci_ub_hi", 
              "eps_zero_lo", "eps_zero_hi")
  n_row <- length(eps_seq)
  n_col <- length(cnames) 
  out <- c(eps_seq, rep(n, n_row), lb, ub, rep(eps_zero$est, n_row),
           ci_lb_lo, ci_ub_hi, rep(eps_zero$ci_lo, n_row), 
           rep(eps_zero$ci_hi, n_row))
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
  sims <- pbreplicate(nsim, sim_fn(n[i], nsplits[i]))
  res[i, "n"] <- n[i]
  
  res[i, "bias_lb"] <- bias(sims[, "lb", ], truth[, "lb"])
  res[i, "bias_ub"] <- bias(sims[, "ub", ], truth[, "ub"])
  res[i, "bias_eps0"] <- bias(sims[1, "eps_zero", ], truth[1, "eps_zero"])
  
  res[i, "rmse_lb"] <- rmse(sims[, "lb", ], truth[, "lb"], n[i])
  res[i, "rmse_ub"] <- rmse(sims[, "ub", ], truth[, "ub"], n[i])
  res[i, "rmse_eps0"] <- rmse(sims[1, "eps_zero", ], truth[1, "eps_zero"], n[i])
  
  res[i, "cvg_reg"] <- coverage(sims[, "ci_lb_lo", ], sims[, "ci_ub_hi", ], 
                                truth[, "lb"], truth[, "ub"])
  res[i, "cvg_eps0"] <- coverage(sims[1, "eps_zero_lo", ], 
                                 sims[1, "eps_zero_hi", ], 
                                 truth[1, "eps_zero"], truth[1, "eps_zero"])

  print(res[1:i, ])
  
  saveRDS(sims, file=paste0("./results/simulation/sims", n[i], ".RData"))
  saveRDS(res, file=paste0("./results/simulation/sim_res", n[i], ".RData"))
}

saveRDS(res, file=paste0("./results/simulation/sim_res.RData"))

print(round(readRDS(paste0("./results/simulation/sim_res.RData")), 2))
