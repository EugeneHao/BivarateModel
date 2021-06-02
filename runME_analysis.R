D = 100 
b = 100
D_each = list(1:25, 26:50, 51:75, 76:100)
library('dplyr')
source("analysis_basefun.R")

# compute MSE and bias ratio: 
result = readRDS("result/simuresult_EBP.rds")

AvMSE_BR(result, pos = 3)   # compute MSE and BR using univariate model
 
AvMSE_BR(result, pos = 4)   # compute MSE and BR using bivariate model 


# compute relative bias:
relative_bias(S = 100, D_each)