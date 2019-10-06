# URL http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.html
# Read dataset into R
# rm(list = ls())
setwd("C:/Users/matte/Dropbox/Causal questions shared with Matteo/Sensitivity Analysis")
library(devtools)
devtools::install("C:/Users/matte/Desktop/sensAteBounds")
library(sensAteBounds)
library(varhandle)
set.seed(1000)
dat <- read.csv("./experiments/data/rhc.csv", header=TRUE)

covariates <- c("age", "sex", "race", "edu", "income",
                
                "cat1", "resp", "card", "neuro", "gastr", "renal", "meta",
                "hema", "seps", "trauma", "ortho",
                
                "ca", "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx",
                "renalhx", "liverhx", "gibledhx", "malighx", "immunhx",
                "transhx", "amihx",
                
                "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1",
                "pafi1", "alb1", "hema1", "bili1", "crea1", "sod1", "pot1", 
                "paco21", "ph1", "wtkilo1", "dnr1",
                
                "ninsclas")

x <- dat[, covariates, drop=FALSE]
# newdf <- x
# for(i in 1:ncol(x)) {
#   if(class(x[, i])=="factor") {
#     new <- to.dummy(x[, i], colnames(x)[i])
#     newdf <- newdf[, !colnames(newdf)%in%colnames(x)[i]]
#     newdf <- cbind(newdf, new[, -1, drop=FALSE])
#   }
# }
x <- newdf
colnames(x) <- paste("x", 1:ncol(x), sep="")
a <- ifelse(dat$swang1 == "No RHC", 0, 1)
y <- ifelse(dat$dth30 == "No", 1, 0)

eps_seq <- seq(0, 0.2, length.out = 200)
delta_seq <- 1
model <- c("x", "xa")
sl.lib <- c("SL.mean", "SL.speedlm", "SL.speedglm", "SL.gam", 
            "SL.ranger", "SL.polymars", "SL.svm")
nuis_fns <- do_crossfit(y=y, a=a, x=x, nsplits=5, outfam=binomial(),
                        treatfam = binomial(), sl.lib=sl.lib)
bounds <- NULL
for(mm in model) {
  res <- get_bound(y=y, a=a, x=x, outfam=NULL, treatfam=NULL, 
                   model=mm, eps=eps_seq, delta=delta_seq, nsplits=NULL, 
                   do_mult_boot=TRUE, do_eps_zero=TRUE,
                   nuis_fns=nuis_fns, alpha=0.05, B=1)
  bounds_next <- res$bounds
  bounds_next <- data.frame(epsilon=res$bounds$eps, delta=res$bounds$delta, 
                            lb=res$bounds$lb, ub=res$bounds$ub, 
                            ci_lo=res$bounds$ci_lo, ci_hi=res$bounds$ci_hi,
                            ci_lo_ptw=res$bounds$ci_lo_ptwise, 
                            ci_hi_ptw=res$bounds$ci_hi_ptwise,
                            ci_lo_ptw_im04=res$bounds$ci_lo_ptwise_im04, 
                            ci_hi_ptw_im04=res$bounds$ci_hi_ptwise_im04,
                            length_bound=res$bounds$ub-res$bounds$lb, 
                            no_zero=res$bound$no_zero,
                            eps_zero=res$bounds$eps_zero,
                            eps_zero_lo=res$bounds$eps_zero_lo,
                            eps_zero_hi=res$bounds$eps_zero_hi)
  bounds_next$model <- mm
  bounds <- rbind(bounds, bounds_next)
}

# sanity checks
stopifnot(sum(round(bounds$length_bound, 10) < 0)==0)
write.csv(bounds, file="./experiments/results/data analysis/rhc_bounds_temp.csv",
row.names=FALSE)


