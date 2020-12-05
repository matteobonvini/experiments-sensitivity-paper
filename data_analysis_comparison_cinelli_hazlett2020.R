# comparison with Cinelli and Hazlett (2020)
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(sensemakr)
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
                "bili1", "alb1", "urin1", "surv2md1",
                
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

# A = 1 means patient underwent RHC
a <- ifelse(dat$swang1 == "No RHC", 0, 1)
# Y = 1 means survival at day 30
y <- ifelse(dat$dth30 == "No", 1, 0)

df <- cbind(data.frame(y=y, a=a), x)

fit <- lm(y ~ ., data=df)

sens_dnr <- sensemakr(fit, treatment = "a", 
                      benchmark_covariates = "dnr1Yes", kd = 2)
sens_miss_adld3p <- sensemakr(fit, treatment = "a", 
                              benchmark_covariates = "is_miss_adld3p", kd = 2)
pdf("./results/data analysis/CinelliHazlett_dnr.pdf")
  plot(sens_dnr, sensitivity.of="t-value", cex.label.text = 1.15,
       cex.axis=1.25)
dev.off()
pdf("./results/data analysis/CinelliHazlett_miss_adld3p.pdf")
plot(sens_miss_adld3p, sensitivity.of="t-value", cex.label.text = 1.15,
     cex.axis=1.5)
dev.off()
robustness_value(fit, covariate="a")
sensitivity_stats(fit, treatment = "a")
ovb_minimal_reporting(sens_dnr)
ovb_minimal_reporting(sens_miss_adld3p)