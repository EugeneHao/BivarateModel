library('rstan')
library('dplyr')
library('parallel')
library('MASS')
library('Rcpp')
library('inline')
library('RcppArmadillo')  

# setwd("F:/Project/SAE/Revision/StrongExample")
source("Bayes_basefun.R")

if(dir.exists("rstan/") == FALSE)
  dir.create("rstan/")


MCsize = 100
D = 100
cores = 50

Pop = 2000
popsize = rep(c(100, 200, 300, 400), each = 25) 
samplesize = popsize * 0.05
fittingsize = popsize - samplesize


result = list()
saveRDS(result, file = "rstan/Bayes_result.rds")


for(i in 1:MCsize)
{
    temp <- Bayes_onesimu(i) 
    result[[i]] <- temp
    saveRDS(temp, file = "rstan/Bayes_result.rds")
}


