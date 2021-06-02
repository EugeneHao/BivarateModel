##### for univariate model: 
multiYprob_uni <- function(mu0) 
{
  prob = (mu0 %*% t(Ychoice)) %>% exp    #n * 2
  rm(mu0)
  sumprob = rowSums(prob)
  prob = sweep(prob, STATS = sumprob, FUN = "/", MARGIN = 1)
  return(prob)
}

## return the log density of Y in univariate case, return a n*m matrix
# mu0:  n * m matrix, Ychoice is 2 * 1 vector 
getgy_uni <- function(mu0, ylabel, m) 
{
  prob = sapply(1:m, FUN = function(x) (mu0[,x] %*% t(Ychoice))) %>% 
    array(., dim = c(dim(mu0)[1], 2, m)) #n * 2 * m
  index = cbind(rep(1:dim(prob)[1], m), rep(ylabel, m), rep(1:m, each = dim(prob)[1]))
  sumprob = apply(exp(prob), MARGIN = c(1,3), FUN = sum) %>% as.vector() %>% log()
  return(matrix(prob[index] - sumprob, ncol = m))
}

# return group log density of Y|X 
# here rand_sample is a D * m matrix 
# update getgy_uni function 
grouplogdensity_uni <- function(ylabel, By_hat, rand_sample, mu_fix, m)
{
  n = sum(samplesize)
  rand_array = rand_sample %>% rep(., rep(samplesize, m)) %>% matrix(., nrow = n)  # n * m matrix
  rm(rand_sample)
  mu = sweep(rand_array, MARGIN = 1, FUN = "+", STATS = mu_fix)  # n * m matrix 
  rm(mu_fix)
  rm(rand_array)
  gy = getgy_uni(mu, ylabel,  m)    # updated 
  rm(mu)
  data.frame(group = rep(1:D, samplesize), gy) -> temp
  result = sapply(1:D, FUN = function(x) temp[temp$group==x, -1] %>% colSums() %>% as.vector) %>% t() 
  return(result)
}

#for multivariate model:
# rand_matrix is D * m matrix 
# update multiYprob_uni
EBPonedraw_uni <-function(data, fitting, By_hat, SigmaV_hat, munew_fix, weight, rand_matrix, m, Z_draw)
{
  onedraw = apply(weight, MARGIN = 1, FUN = function(x) sample(1:m, size = 1, prob = x))
  rm(weight)
  rand_one = sapply(1:D,  FUN = function(x) rand_matrix[x, onedraw[x]]) # scalar 
  rm(rand_matrix)
  rm(onedraw)
  mu_new = munew_fix + rep(rand_one, predictsize) # (N-n) * 1
  rm(munew_fix)
  Yprob_new = multiYprob_uni(mu_new)  #updated 
  Y_draw = sapply(1:dim(Yprob_new)[1], FUN = function(x) sample(1:2, size = 1, prob = Yprob_new[x,]))
  n = dim(Yprob_new)[1]
  Yid = Y_draw %>% as.factor()
  
  Y_draw = Ychoice[Y_draw]
  rm(mu_new)
  rm(Yprob_new)
  data.frame(fitting, Y = Y_draw, Z = Z_draw, Yid = Yid) %>% 
    rbind(., data.frame(data[,1:(num_Zcoef +4)])) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z),
                                  y10mean = mean(1-Y), y11mean = mean(Y), 
                                  z1y10 = sum((1-Y)*Z)/sum(1-Y), z1y11 = sum(Y*Z)/sum(Y)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


## get the EBP entire results, a D * P * t array
#logdens is a m * D  matrix
# we need to update grouplogdensity_uni, EBPonedraw_uni
EBPfun_uni <- function(data, fitting, By_hat, SigmaV_hat, m, tt, Zdraw)
{
  aux = cbind(1, data$X2)
  mu_fix = aux %*% By_hat   # n * 1
  Y = data$Y
  ylabel = data$Yid
  munew_fix = cbind(1,fitting$X2) %*% By_hat
  rand_matrix= rnorm(m*D, sd = sqrt(SigmaV_hat)) %>% matrix(., nrow = D)  # D * m * 1 
  logdens = grouplogdensity_uni(ylabel, By_hat, rand_matrix,  mu_fix, m)   # updated 
  dens = (logdens - median(logdens)) %>% round(3) %>% exp()
  weight = sweep(dens, STATS= rowSums(dens), FUN = "/" , MARGIN = 1)
  
  temp = sapply(1:tt, FUN = function(x) EBPonedraw_uni(data, fitting, By_hat, SigmaV_hat, 
                                                       munew_fix, weight, rand_matrix, m, Zdraw)) %>% 
    array(., dim = c(D, 5, tt))
  
  return(temp)
}
