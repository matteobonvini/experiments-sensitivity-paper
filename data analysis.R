############################################################
## Analyze data of Connors et al (1996) as in Section 4.2 ##
############################################################
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(varhandle)
library(devtools)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # if warnings are causing erros and they are not important one may use this
# devtools::install_github("matteobonvini/sensitivitypuc", force = TRUE)
library(sensitivitypuc)

options(stringsAsFactors = TRUE)

set.seed(1000)

data_url <- "http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv"
dat <- read.csv(data_url, header = TRUE)

covariates <- c("age", "sex", "race", "edu", "income", 
                "ninsclas", "cat1", "cat2", 
                "resp", "card", "neuro", "gastr", "renal", "meta",
                "hema", "seps", "trauma", "ortho",
                
                "adld3p", "das2d3pc", "dnr1", "ca",  "aps1", "scoma1", 
                "wtkilo1", "temp1", "meanbp1", "resp1", "hrt1", "pafi1", 
                "paco21", "ph1", "wblc1", "hema1", "sod1", "pot1", "crea1", 
                "bili1", "alb1", "urin1",
                
                "cardiohx", "chfhx", "dementhx", 
                "psychhx", "chrpulhx", "renalhx", "liverhx", "gibledhx", 
                "malighx", "immunhx", "transhx", "amihx")

x <- dat[, covariates]
miss_covs <- covariates[apply(x, 2, function(x) sum(is.na(x)) > 0)]

missdat <- as.data.frame(matrix(0, ncol = length(miss_covs), nrow = nrow(x),
                         dimnames = list(NULL, paste0("is_miss_", miss_covs))))

for(varname in miss_covs) {
  idx <- which(is.na(x[, varname]))
  if(class(x[, varname]) == "factor") {
    levels(x[, varname]) <- c(levels(x[, varname]), "missing")
    x[idx, varname] <- "missing"
    missdat[, paste0("is_miss_", varname)] <- NULL
  } else {
    x[idx, varname] <- 0
    missdat[idx, paste0("is_miss_", varname)] <- 1
  }
}

x <- cbind(x, missdat)
# Fix factors or otherwise SuperLearner complains
newdf <- x
for(i in 1:ncol(x)) {
  if(class(x[, i]) == "factor") {
    tmp <- to.dummy(x[, i], colnames(x)[i])
    newdf <- newdf[, !colnames(newdf) %in% colnames(x)[i]]
    newdf <- cbind(newdf, tmp[, -1, drop = FALSE])
  }
}
x <- newdf
rm(newdf)
rm(tmp)

# Avoid SuperLearner conflicts due to variable names
colnames(x) <- paste("x", 1:ncol(x), sep = "")
# A = 1 means patient underwent RHC
a <- ifelse(dat$swang1 == "No RHC", 0, 1)
# Y = 1 means survival at day 30
y <- ifelse(dat$dth30 == "No", 1, 0)

# unadjusted OR (using morality as outcome)
mean(y[a == 0]) / (1 - mean(y[a == 0])) * (1 - mean(y[a == 1])) / mean(y[a == 1]) 
# unadjusted RD
mean(y[a == 0]) - mean(y[a == 1])

# Select values for prop of unmeasured confounding at which evaluate bounds
eps_seq <- seq(0, 0.5, 0.0001)
delta_seq <- c(0.25, 0.50, 0.75, 1)
# Select model, "x" = S \ind (Y, A) | X, "xa" = S \ind Y | (X, A) 
model <- c("x", "xa")

# Select confidence level and # rademachers for multiplier bootstrap
alpha <- 0.05
B <- 10000
# Select SuperLearner Library
sl.lib <- c("SL.mean", "SL.speedlm", "SL.speedglm", "SL.gam", "SL.svm", 
            "SL.polymars", "SL.ranger")
nsplits <- 5
# Estimate Regression functions once for both model "x" and model "xa"
# There is a "non-list contrasts argument ignored" warning from SL gam library
# coming from model.matrix, which I think can be ignored. 
# system.time({
#   nuis_fns <- do_crossfit(y = y, a = a, x = x, ymin = 0, ymax = 1,
#                           nsplits = nsplits, outfam = binomial(),
#                           treatfam = binomial(), sl.lib = sl.lib,
#                           do_parallel = TRUE, ncluster = 3)
# })
rch_file_name <- paste0("./results/data analysis/nuis_fns_rhc_", nsplits,
                        "fold.RData")
# saveRDS(nuis_fns, file = rch_file_name)
nuis_fns <- readRDS(rch_file_name)

system.time({
res_x <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, model = "x", 
                   eps = eps_seq, delta = delta_seq, nuis_fns = nuis_fns, 
                   alpha = alpha, B = B, do_mult_boot = TRUE,
                   do_eps_zero = TRUE, do_rearrange = TRUE)
})

system.time({
res_xa <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, model = "xa", 
                    eps = eps_seq, delta = delta_seq, nuis_fns = nuis_fns, 
                    alpha = alpha, B = B, do_mult_boot = TRUE,
                    do_eps_zero = TRUE, do_rearrange = TRUE)
})

# saveRDS(res_x, file = paste0("./results/data analysis/res_x_", nsplits,
                             # "fold.RData"))
# saveRDS(res_xa, file = paste0("./results/data analysis/res_xa_", nsplits,
                             # "fold.RData"))

# res_x <- readRDS(paste0("./results/data analysis/res_x_", nsplits, "fold.RData"))
# res_xa <- readRDS(paste0("./results/data analysis/res_xa_", nsplits, "fold.RData"))

bound_x <- res_x$bounds
bound_xa <- res_xa$bounds
eps_zero_x <- res_x$eps_zero
eps_zero_xa <- res_xa$eps_zero   
round(100 * bound_x[1, , "1"], 2)
round(100 * bound_xa[1, , "1"], 2)
round(100 * eps_zero_x, 2)
round(100 * eps_zero_xa, 2)

# Export result in a user-friendly format
out <- NULL
for(i in 1:length(delta_seq)) {
  # Loop through the values of deltas and append
  dd <- delta_seq[i]
  new_bound_x <- cbind(eps_seq, dd, bound_x[, , as.character(dd)], 
                   eps_zero_x$est[i], eps_zero_x$ci_lo[i], eps_zero_x$ci_hi[i], 
                   "x")
  new_bound_xa <- cbind(eps_seq, dd, bound_xa[, , as.character(dd)], 
                    eps_zero_xa$est[i], eps_zero_xa$ci_lo[i], 
                    eps_zero_xa$ci_hi[i], "xa")
  bound <- rbind(new_bound_x, new_bound_xa)
  out <- rbind(out, bound)
}
colnames(out) <- c("epsilon", "delta", colnames(res_x$bounds[, , 1]), 
                   "eps_zero", "eps_zero_lo", "eps_zero_hi", "model")
write.csv(out, "./results/data analysis/rhc_bounds.csv", row.names = FALSE)
