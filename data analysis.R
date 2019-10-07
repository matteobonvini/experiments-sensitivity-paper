##########################################################
## Analyze data of Connors et al (1996) as in Section 5 ##
##########################################################
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(devtools)
devtools::install("C:/Users/matte/Desktop/sensAteBounds")
library(sensAteBounds)
library(varhandle)

set.seed(1000)
# Reference: http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.html
dat <- read.csv("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.csv", 
                header = TRUE)

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
print(miss_covs)
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
eps_seq <- seq(0, 0.2, 0.001)
# Select model, "x" = S \ind (Y, A) | X, "xa" = S \ind Y | (X, A) 
model <- c("x", "xa")

# Select SuperLearner Library
sl.lib <- c("SL.mean", "SL.speedlm", "SL.speedglm", "SL.gam", 
            "SL.ranger", "SL.polymars", "SL.svm")
# Estimate Regression functions once for both model "x" and model "xa"
# There is a "non-list contrasts argument ignored" warning from SL gam library
# coming from model.matrix, which I think can be ignored. 
nuis_fns <- do_crossfit(y = y, a = a, x = x, nsplits = 5, outfam = binomial(),
                        treatfam = binomial(), sl.lib = sl.lib)
saveRDS(nuis_fns, file = "./results/nuis_fns_rhc.RData")
nuis_fns <- readRDS("./results/nuis_fns_rhc.RData")
bounds <- NULL
for(mm in model) {
  # Loop thru models and evaluate bounds curves
  res <- get_bound(y=y, a=a, x=x, outfam=NULL, treatfam=NULL, 
                   model=mm, eps=eps_seq, delta=1, nsplits=NULL, 
                   do_mult_boot=TRUE, do_eps_zero=TRUE,
                   nuis_fns=nuis_fns, alpha=0.05, B=10000)
  bounds_next <- data.frame(epsilon=res$bounds$eps, delta=res$bounds$delta, 
                            lb=res$bounds$lb, ub=res$bounds$ub, 
                            ci_lo=res$bounds$ci_lo, ci_hi=res$bounds$ci_hi,
                            ci_lo_ptw=res$bounds$ci_lo_ptwise, 
                            ci_hi_ptw=res$bounds$ci_hi_ptwise,
                            ci_lo_ptw_im04=res$bounds$ci_lo_ptwise_im04, 
                            ci_hi_ptw_im04=res$bounds$ci_hi_ptwise_im04,
                            length_bound=res$bounds$ub-res$bounds$lb, 
                            no_zero=res$bound$no_zero,
                            eps_zero=res$eps_zero$est,
                            eps_zero_lo=res$eps_zero$ci_lo,
                            eps_zero_hi=res$eps_zero$ci_hi)
  bounds_next$model <- mm
  bounds <- rbind(bounds, bounds_next)
}

write.csv(bounds, "./results/data analysis/rhc_bounds.csv", row.names = FALSE)


