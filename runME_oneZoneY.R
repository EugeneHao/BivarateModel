#This code is used for one continuous variable and one categorical variable: 

library('dplyr')
library('tidyr')
library('mvtnorm')
library('lme4')
library('parallel')
library('MASS')
library('Rcpp')
library('inline')
library('RcppArmadillo')  

# you might need to set your own directory 
setwd("/work/STAT/hao123/GauMultJoint/revision2")
source("basefun.R")                            #source three codes 
source("basefun_uni.R")
source("cpp.R")

if(dir.exists("result/") == FALSE)
  dir.create("result/")
if(dir.exists("data/") == FALSE)
  dir.create("data/")


##############################################
# important settings: 
generate = T                      # if TRUE, regenerate the popolation and sample data
noninfor = T                      # if FALSE, use informative sampling 
doestimate = T                    # if TRUE, estimate model parameters, if FALSE, load the saved estimates 
update_all = F                    # if TRUE, use full EM algorithm to update coefficient estimates
doEBP = T                         # if TRUE, compute EBP
usedev = F                        # if TRUE, use derivative function for optimazation 
wronguse = F                      # if TRUE, use informative sampling data, but ignore the weights 
bootstrap = F                     # if TRUE, use bootstrap to estimate MSE 
useuni = T                        # if TRUE, also use univariate independence model 
case = 1                          # example for Moderate (+) in Section 4.1 


# create your file to save the result: 
if(doEBP == T)
  filename = paste("result/simuresult_EBP", ".rds", sep = "")   
if(doEBP == F & doestimate == T)
  filename = paste("result/simuresult_estimate", ".rds", sep = "")   
if(doEBP == F & doestimate == F)
  stop("please set doEBP = T or doestimate = T")


result = list()
saveRDS(result, file = filename)


##############################################
# computation settings: 
m = 200               # Monte Carlo size
t = 200               # EBP size (equal to n_eff)
b = 100               # bootstrap size
MCsize = 100          # Monte Carlo size
seed = 1:MCsize       # set the seed

EMsize = 20           # iteraton times for EM alrorithm
cores = 16            # the number of cores used 

################################################
# model settings: 
D = 100                                            # group size 
Pop = rep(c(100, 200, 300, 400), each = D/4)       # population size 
samplesize = Pop * 0.05                            # samplesize
popsize = Pop                                     
fittingsize = popsize - samplesize                  
predictsize = popsize - samplesize 

# model paramters:           
if(case == 1)   # Moderate (+)
{
  Bz = c(0, 1)                       # beta_z
  By = c(-2, 0.5)                    # beta_y
  P = 0.5                            # rho_y
  SigmaZ = 1                         # sigma^2
  rho_uv = 0.4
}       

if(case == 2)   # Moderate (-)
{
  Bz = c(0, 1)                       # beta_z
  By = c(-2, 0.5)                    # beta_y
  P = 0.5                            # rho_y
  SigmaZ = 1                         # sigma^2
  rho_uv = -0.4
}       

if(case == 3)   # Strong (rho)
{
  Bz = c(0, 1)                       # beta_z
  By = c(-4, 0.5)                    # beta_y
  P = 2                              # rho_y
  SigmaZ = 0.75                      # sigma^2
  rho_uv = 0.4
}       

if(case == 4)   # Strong (x)
{
  Bz = c(-1, 1, 1, 1.5)              # beta_z
  By = c(-4, 0.5)                    # beta_y
  P = 2                              # rho_y
  SigmaZ = 0.75                      # sigma^2
  rho_uv = 0.4
}       


num_Zcoef = length(Bz) - 1                 
num_Ycoef = length(By) - 1          
fz = paste("Z~", paste(paste("X", 2:(num_Zcoef+1), sep = ""), collapse = "+"), "+ Yid + (1|group)", sep = "")
fz_uni = paste("Z~", paste(paste("X", 2:(num_Zcoef+1), sep = ""), collapse = "+"), " + (1|group)", sep = "")


SigmaUV = diag(c(1, 0.5)) %*% 
  matrix(c(1, rho_uv,
           rho_uv, 1), nrow = 2) %*%
  diag(c(1, 0.5))               #variance of random effects
Ychoice = c(0,1)

# you can do the parallel computation or use the for loop
result = list()
for(x in 1: MCsize)
{
   temp <-onestep(m, t, b, seed[x], EMsize)
   result[[x]] <- temp
   saveRDS(result, file = filename)
}

# mclapply(seed, FUN = function(x) onestep(m, t, b, seed[x], EMsize), mc.cores = cores) -> result
# saveRDS(result, file = filename)

