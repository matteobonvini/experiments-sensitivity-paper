##########################################################
## Analyze data of Connors et al (1996) as in Section 5 ##
##########################################################
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(devtools)
library(varhandle)
library(sensitivitypuc)

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

# Exclude covariates with missing values
x <- dat[, covariates]
miss_covs <- covariates[apply(x, 2, function(x) sum(is.na(x)) > 0)]
print(paste("Discarding covariates:", paste(miss_covs, collapse = ", ")))
x <- x[, !covariates %in% miss_covs]

# Fix factors or otherwise SuperLearner complains
newdf <- x
for(i in 1:ncol(x)) {
  if(class(x[, i])=="factor") {
    new <- to.dummy(x[, i], colnames(x)[i])
    newdf <- newdf[, !colnames(newdf)%in%colnames(x)[i]]
    newdf <- cbind(newdf, new[, -1, drop=FALSE])
  }
}
x <- newdf
rm(newdf)
rm(new)

# Avoid SuperLearner conflicts due to variable names
colnames(x) <- paste("x", 1:ncol(x), sep = "")
# A = 1 means patient underwent RHC
a <- ifelse(dat$swang1 == "No RHC", 0, 1)
# Y = 1 means survival at day 30
y <- ifelse(dat$dth30 == "No", 1, 0)

# Select values for prop of unmeasured confounding at which evaluate bounds
eps_seq <- seq(0, 0.6, 0.0001)
delta_seq <- c(0.25, 0.50, 0.75, 1)
# Select model, "x" = S \ind (Y, A) | X, "xa" = S \ind Y | (X, A) 
model <- c("x", "xa")

# Select confidence level and # rademachers for multiplier bootstrap
alpha <- 0.05
B <- 10000
# Select SuperLearner Library
sl.lib <- c("SL.mean", "SL.speedlm", "SL.speedglm", "SL.gam",
            "SL.ranger", "SL.polymars", "SL.svm")
# Estimate Regression functions once for both model "x" and model "xa"
# There is a "non-list contrasts argument ignored" warning from SL gam library
# coming from model.matrix, which I think can be ignored. 
nuis_fns <- do_crossfit(y = y, a = a, x = x, ymin = 0, ymax = 1, nsplits = 5,
                        outfam = binomial(), treatfam = binomial(),
                        sl.lib = sl.lib, do_parallel = TRUE, ncluster = 3,
                        show_progress = FALSE)
saveRDS(nuis_fns, file = "./results/data analysis/nuis_fns_rhc.RData")
nuis_fns <- readRDS("./results/data analysis/nuis_fns_rhc.RData")

res_x <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, model = "x", 
                   eps = eps_seq, delta = delta_seq, nuis_fns = nuis_fns, 
                   alpha = alpha, B = B, do_mult_boot = TRUE,
                   do_eps_zero = TRUE, do_rearrange = TRUE)

res_xa <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, model = "xa", 
                    eps = eps_seq, delta = delta_seq, nuis_fns = nuis_fns, 
                    alpha = alpha, B = B, do_mult_boot = TRUE,
                    do_eps_zero = TRUE, do_rearrange = TRUE)

bound_x <- res_x$bounds
bound_xa <- res_xa$bounds
eps_zero_x <- res_x$eps_zero
eps_zero_xa <- res_xa$eps_zero   

# Export result in a user-friendly format
out <- NULL
for(i in 1:length(delta_seq)) {
  # Loop through the values of deltas and append
  dd <- delta_seq[i]
  bound_x <- cbind(eps_seq, dd, res_x$bounds[, , as.character(dd)], 
                   eps_zero_x$est[i], eps_zero_x$ci_lo[i], eps_zero_x$ci_hi[i], 
                   "x")
  bound_xa <- cbind(eps_seq, dd, res_xa$bounds[, , as.character(dd)], 
                    eps_zero_xa$est[i], eps_zero_xa$ci_lo[i], 
                    eps_zero_xa$ci_hi[i], "xa")
  bound <- rbind(bound_x, bound_xa)
  out <- rbind(out, bound)
}
colnames(out) <- c("epsilon", "delta", colnames(res_x$bounds[, , 1]), 
                   "eps_zero", "eps_zero_lo", "eps_zero_hi", "model")
write.csv(out, "./results/data analysis/rhc_bounds.csv", row.names = FALSE)