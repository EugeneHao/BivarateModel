# functions 
multiYprob <- function(P0, SigmaZ0, delta0, mu0) 
{
   gx = as.vector(P0 * delta0 + mu0 + SigmaZ0 * P0^2/2)
   prob = cbind(1/(1+exp(gx)), exp(gx)/(1+exp(gx)))   # 0 first, then 1
  return(prob)
}

# update rcppYZdraw
rcppEBPonedraw <-function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, deltanew_fix, munew_fix, weight, rand_matrix, m)
{
  onedraw = apply(weight, MARGIN = 1, FUN = function(x) sample(1:m, size = 1, prob = x))  # weight D * m
  rm(weight)
  rand_one = sapply(1:D,  FUN = function(x) rand_matrix[x, onedraw[x], ]) %>% t()
  rm(rand_matrix)
  rm(onedraw)
  predictsize = table(fitting$group) %>% as.vector() 
  data_unique = data %>% distinct()

  delta_new = deltanew_fix + rep(rand_one[,1], predictsize)
  mu_new = munew_fix + rep(rand_one[,2], predictsize)
  rm(deltanew_fix)
  rm(munew_fix)
  Yprob_new = multiYprob(P_hat, SigmaZ_hat, delta_new, mu_new)  # (N-n) * 2

     Y = lapply(1:dim(Yprob_new)[1], FUN = function(x) sample(1:2, size = 1, prob = Yprob_new[x,])) %>% unlist()
     Yid = Y  %>% as.factor()
     Y_draw = Ychoice[Y]  
     condZmean = delta_new + SigmaZ_hat * Y_draw * P_hat 
     Z_draw = rnorm(sum(predictsize), mean = condZmean, sd = sqrt(SigmaZ_hat)) 

  data.frame(fitting, Y = Y_draw, Z = Z_draw, Yid =Yid) %>% 
    rbind(., data.frame(data_unique[,1:(num_Zcoef + 4)])) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z, na.rm = T), y10mean = mean(1-Y,na.rm = T), y11mean = mean(Y, na.rm = T), 
                                  z1y10 = sum((1-Y)*Z, na.rm = T)/sum(1-Y, na.rm = T), z1y11 = sum(Y*Z, na.rm = T)/sum(Y, na.rm = T)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}

rcppEBPfun <- function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat,  SigmaUV_hat, m, tt)
{
  aux = cbind(1, data[,2:(1+num_Zcoef)]) %>% as.matrix()
  if(noninfor == FALSE)
  {
    lw = lm(log(weights) ~ 0 + Z + Y + as.factor(data$group), data = data)
    xi = lw$coef
    delta_fix = aux %*% (Bz_hat + c(xi[1]*SigmaZ_hat, 0))
    mu_fix = aux[,1:2] %*% (By_hat + c(xi[2], 0))
  }
  if(noninfor == TRUE)
  {
    delta_fix = aux %*% Bz_hat
    mu_fix = aux[,1:2] %*% By_hat
  }
  deltanew_fix = as.matrix(cbind(1,fitting[, 2:(1+num_Zcoef)])) %*% Bz_hat
  munew_fix = cbind(1,fitting$X2) %*% By_hat

  rand_matrix= rmvnorm(m*D, sigma = SigmaUV_hat) %>% array(., dim = c(D, m, 2))
  
  # logdens: m * D
  logdens = rcppgroupdens(m0 = m, D0 = D, P0 = P_hat, sigmaZ0 = SigmaZ_hat, Z0 = data$Z, Y0 = data$Y, 
                          rand_matrix = rand_matrix, delta_fix0 = delta_fix, mu_fix0 = mu_fix, id0 = data$group, subid0 = data$subgroupid)
  dens = (logdens - median(logdens)) %>% exp() 
  weight = sweep(dens, STATS= colSums(dens), FUN = "/" , MARGIN = 2) %>% t()
  
  temp = sapply(1:tt, FUN = function(x) rcppEBPonedraw(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, 
                                                       deltanew_fix, munew_fix, weight, rand_matrix, m)) %>% 
    array(., dim = c(D, 5, tt))
  
  return(temp)
}


## Simulate sample data for both true data and bootstrap
simuldata <- function(auxillary, Bz_hat, By_hat,  P_hat, SigmaZ_hat, SigmaUV_hat, groupsize)
{
  randef = rmvnorm(D, sigma = SigmaUV_hat) %>%
    rep(., rep(groupsize,2)) %>% matrix(ncol = 2) %>% as.data.frame() %>% 
    "names<-"(c("u1", "v1"))   # sum(groupsize) * 2
  X = cbind(1, auxillary)
  delta_draw = (X %*% Bz_hat + randef[, 1]) %>% as.matrix()
  mu_draw = (X[,1:2] %*% By_hat + randef[, 2]) %>% as.matrix()
  group = rep(c(1:D), groupsize)   
  Yprob = multiYprob(P_hat, SigmaZ_hat, delta_draw, mu_draw)  # n * 2
  Y = lapply(1:dim(Yprob)[1], FUN = function(x) sample(1:2, size = 1, prob = Yprob[x,])) %>% unlist()
  Yid = Y  %>% as.factor()
  Y = Ychoice[Y]  # %>% 'colnames<-'(paste("Y",  1:2, sep = ""))
  condZmean = delta_draw + SigmaZ_hat * Y * P_hat 
  Z = rnorm(sum(groupsize), mean = condZmean, sd = sqrt(SigmaZ_hat)) 
  popdata = data.frame(group, auxillary, Y, Z, Yid) %>% 
	"names<-"(c('group', paste("X", 2:(1+num_Zcoef), sep=""), "Y", "Z", "Yid"))
  return(popdata)
}

# one replicate for bootstrap method
oneboot <- function(data, fitting, Bz0, By0, SigmaZ0, P0, SigmaUV0, m, tt, EMsize)
{
  popaux <- rbind(data[,1:(1+num_Zcoef)], fitting) %>% arrange(., group)
  X0 = popaux[,2:(1+num_Zcoef)]
  truedata <- simuldata(X0, Bz0, By0, P0, SigmaZ0, SigmaUV0, popsize) 
  
  if(noninfor == FALSE)
  {
    xi = rnorm(sum(popsize))
    xi[xi > 2] = 2
    xi[xi <-2] = -2
    numeriator = exp(-truedata$Z/3 + truedata$Y/2 + xi/5)
    denomiator = cbind(truedata, numeriator) %>% group_by(group) %>% summarize(sumnum = sum(numeriator)) %>% "$"(sumnum) %>% rep(., popsize)
    piij = rep(samplesize, popsize) * numeriator/denomiator
    truedata$weights = piij
  }
  sampledata <- left_join(data[,1:(1+num_Zcoef)], truedata) # include weights under informative sampling 
  sampledata$subgroupid = rep(1, dim(sampledata)[1])
  
  # model parameter estimate 
  lm_z = lmer(as.formula(fz), data = sampledata)
  Bz_hat = fixef(lm_z)[1:1:(1+num_Zcoef)]
  U_hat = ranef(lm_z)$group
  SigmaZ_hat = VarCorr(lm_z) %>% as.data.frame() %>% '['(2,4)
  P_hat1 = fixef(lm_z)[2+num_Zcoef]/SigmaZ_hat
  prec_P1 = 1/(vcov(lm_z)[2+num_Zcoef, 2+num_Zcoef]/SigmaZ_hat^2)
  
  glm_y = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata, family = binomial, nAGQ = 20)
  By_hat = glm_y %>% fixef %>% "["(1:2)
  
  V_hat <- rep(0, 40)
  V_hat[rownames(ranef(glm_y)$group ) %>% as.numeric()] <- ranef(glm_y)$group %>% unlist
  
  rand_hat = cbind(U_hat, V_hat)
  smallsd = sqrt(c(0.0001, 0.0008))/2
  if(isSingular(lm_z) == TRUE)
    rand_hat[,1] = rnorm(D, sd = smallsd[1])
  if(isSingular(glm_y) == TRUE)
    rand_hat[,2] = rnorm(D, sd = smallsd[2])
  SigmaUV_hat = rand_hat %>% cov
  P_hat2 = fixef(glm_y)[3]
  prec_P2 = 1/diag(vcov(glm_y))[3]
  
  P_hat = 1/(prec_P1 + prec_P2) * (prec_P1 * P_hat1 + prec_P2 * P_hat2)
  
  if(noninfor == FALSE && wronguse == FALSE)
  {
    lw = lm(log(weights) ~ 0 + Z + Y + as.factor(sampledata$group), data = sampledata)
    zeta = lw$coef
    By_hat[1] = By_hat[1] - zeta[2]
    Bz_hat[1] = Bz_hat[1] - zeta[1]*SigmaZ_hat
  }
  
  # MCEM 
  initialest = c(Bz_hat, By_hat, P_hat, SigmaZ_hat, as.vector(SigmaUV_hat))
  for(k in 1:200)
    SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m = 2000)
  
  for(j in 1:EMsize)
  {
    for(k in 1:5)
      SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m = 2000)
    if(update_all == TRUE)
    {
      rcppupdatefix(sampledata, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, maxiter = 50, m = 3000) -> temp
      Bz_hat = temp[1:(1+num_Zcoef)] 
      By_hat = temp[(2+num_Zcoef):(2+num_Zcoef+num_Ycoef)]
      P_hat = temp[3+num_Zcoef+num_Ycoef]
      SigmaZ_hat = temp[4+num_Zcoef+num_Ycoef]
    }
  }
  
  # bootstrap mean and variance 
  rcppEBPfun(data, fitting, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, m, tt) -> temp  # D * 5 * tt
  bootmean = apply(temp, MARGIN = c(1,2), FUN = mean) %>% as.vector()
  bootvar = apply(temp, MARGIN = c(1,2), FUN = var) %>% as.vector()
  boottrue = truedata %>% group_by(group) %>% 
              summarise(z1mean = mean(Z, na.rm = T), y10mean = mean(1-Y,na.rm = T), y11mean = mean(Y, na.rm = T), 
                        z1y10 = sum((1-Y)*Z, na.rm = T)/sum(1-Y, na.rm = T), 
		    z1y11 = sum(Y*Z, na.rm = T)/sum(Y, na.rm = T)) 
  return(cbind(bootmean, bootvar, boottrue[ ,-1] %>% unlist() %>% as.vector()))
}



# Main Function
onestep <- function(m, t, b, x, EMsize)
{
  if(generate == TRUE)   
  {
    if(num_Zcoef == 3)
    {
      X2 = runif(sum(popsize), min = 1.5, max  = 3.5)    
      X3 = runif(sum(popsize), min = 1.5, max  = 3.5)    
      X4 = (X3 - mean(X3))^2 - mean((X3 - mean(X3))^2)
      X = cbind(X2, X3, X4)
    }
    if(num_Zcoef == 1)
    {
      X2 = runif(sum(popsize), min = 1.5, max  = 3.5)   
      X = X2
    } 
    truedata = simuldata(X, Bz, By, P, SigmaZ, SigmaUV, popsize)   #   group, X, Y, Z, Yid
    
    saveRDS(truedata, file = paste("data/SAE_pop_", x, ".rds", sep = "")) 
    
    if(noninfor == TRUE)
    {
      # select sample data using SRS
      mapply(sample_n, split(truedata, truedata$group), samplesize, SIMPLIFY = FALSE) %>%
        bind_rows() -> sampledata
      sampledata$subgroupid = rep(1, dim(sampledata)[1])   # 1 denote observing both Y and Z, we do not consider missing data 
      # save simulation data 
      saveRDS(sampledata, file = paste("data/SAE_sample_", x, ".rds", sep = ""))  #   group, X2, X3, X4, Y, Z, Yid, subgroupid
    }
    if(noninfor == FALSE)
    {
       # construct sample weights 
       xi = rnorm(sum(popsize))
       xi[xi > 2] = 2
       xi[xi <-2] = -2
       numeriator = exp(-truedata$Z/3 + truedata$Y/2 + xi/5)
       denomiator = cbind(truedata, numeriator) %>% group_by(group) %>% summarize(sumnum = sum(numeriator)) %>% "$"(sumnum) %>% rep(., popsize)
       piij = rep(samplesize, popsize) * numeriator/denomiator
       truedata$weights = piij
       splitdata = split(truedata, truedata$group)
       splitweight = split(piij, truedata$group)
       # select sample data using informative sampling 
       sampledata = splitdata[[1]][sample(1:popsize[1], samplesize[1], replace = TRUE, prob = splitweight[[1]]), ]
       for (i in 2:D) 
          sampledata = rbind(sampledata, 
              splitdata[[i]][sample(1:popsize[i], samplesize[i], replace = TRUE, prob = splitweight[[i]]), ])
       sampledata$subgroupid = rep(1, dim(sampledata)[1])
       # save simulation data 
       saveRDS(sampledata, file = paste("data/SAE_sample_infor_", x, ".rds", sep = "")) 	
    }
  }

  if(generate == FALSE)
  {
     if(noninfor == TRUE)
     { 
       sampledata = readRDS(file = paste("data/SAE_sample_", x, ".rds", sep = ""))  
     }
     if(noninfor == FALSE || wronguse == TRUE)   
     {
       sampledata = readRDS(file = paste("data/SAE_sample_infor_", x, ".rds", sep = ""))  
     }
    truedata = readRDS(file = paste("data/SAE_pop_", x, ".rds", sep = ""))  
  }


if(doestimate == TRUE)
{
  #first we use Z|Y
  lm_z = lmer(as.formula(fz), data = sampledata)
  Bz_hat = fixef(lm_z)[1:(1+num_Zcoef)]  
  U_hat = ranef(lm_z)$group
  SigmaZ_hat = VarCorr(lm_z) %>% as.data.frame() %>% '['(2,4)
  P_hat1 = as.vector(fixef(lm_z)[2+num_Zcoef])/SigmaZ_hat
  prec_P1 = solve(vcov(lm_z)[2+num_Zcoef, 2+num_Zcoef]/SigmaZ_hat^2)
  
  #next we use Y_i|Z, Y_{-i}  
  glm_y = glmer(Yid ~ X2 + Z  + (1|group), data = sampledata, family = binomial, nAGQ = 20)
  By_hat = glm_y %>% fixef %>% "["(c(1:2))

  V_hat <- rep(0, D)
  V_hat[rownames(ranef(glm_y)$group ) %>% as.numeric()] <- ranef(glm_y)$group %>% unlist
  
  rand_hat = cbind(U_hat, V_hat)
  # use small values as estimates of random effects if the model is singular
  smallsd = sqrt(c(0.0001, 0.0008))/2
  if(isSingular(lm_z) == TRUE)
    rand_hat[,1] = rnorm(D, sd = smallsd[1])
  if(isSingular(glm_y) == TRUE)
    rand_hat[,2] = rnorm(D, sd = smallsd[2])

  SigmaUV_hat = rand_hat %>% cov
  
  P_hat2 = fixef(glm_y)[3]
  prec_P2 = diag(vcov(glm_y))[3]^-1
  
  P_hat = as.vector(1/(prec_P1 + prec_P2) * (prec_P1 * P_hat1 + prec_P2 * P_hat2))
  
  # Monte Carlo EM: 
  initialest = c(Bz_hat, By_hat, P_hat, SigmaZ_hat, as.vector(SigmaUV_hat))
  for(k in 1:200)
    SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m = 2000)
    
  for(j in 1:EMsize)
  {
    for(k in 1:5)
    SigmaUV_hat = rcppSigmaUVfun(sampledata, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m = 2000)
    if(update_all == TRUE)
    {
      rcppupdatefix(sampledata, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, maxiter = 100, m = 3000) -> temp
      Bz_hat = temp[1:(1+num_Zcoef)] 
      By_hat = temp[(2+num_Zcoef):(2+num_Zcoef+num_Ycoef)]
      P_hat = temp[3+num_Zcoef+num_Ycoef]
      SigmaZ_hat = temp[4+num_Zcoef+num_Ycoef]
    }
  }
  
  # under informative case, adjust beta_z0, beta_y0 and rho if needed 
  if(noninfor == FALSE && wronguse == FALSE)
  {
    lw = lm(log(weights) ~ 0 + Z + Y + as.factor(sampledata$group), data = sampledata)
    zeta = lw$coef
    By_hat[1] = By_hat[1] - zeta[2]
    Bz_hat[1] = Bz_hat[1] - zeta[1]*SigmaZ_hat
  }
}

# compute empirical Bayes predictors and use bootstrap methods to estimate MSE: 
if(doEBP == TRUE)
{
  if (doestimate == FALSE)    # if we do not do estimation, load model parameter estimates first
  {
     paramest <- readRDS("result/simuresult_estimate.rds")[[x]]   # check file name!
     Bz_hat <- paramest[2, 1:(1+num_Zcoef)] 
     By_hat <- paramest[2, (2+num_Zcoef):(2+num_Zcoef+num_Ycoef)]
     P_hat <- paramest[2, 3+num_Zcoef+num_Ycoef]
     SigmaZ_hat <- paramest[2, 4+num_Zcoef+num_Ycoef]
     SigmaUV_hat <- paramest[2,(5+num_Zcoef+num_Ycoef):(8+num_Zcoef+num_Ycoef)] %>% matrix(., nrow = 2)
  }

  unsampledata = anti_join(truedata, sampledata) %>% "["(,1:(1+num_Zcoef))
  
  #true mean 
  truedata %>% group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y), y11mean = mean(Y), 
                                             z1y10 = sum((1-Y)*Z)/sum(1-Y), z1y11 = sum(Y*Z)/sum(Y)) -> truemean
  
  #use direct method to compute small area means and domain means 
  sampledata %>% group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y), y11mean = mean(Y), 
                                               z1y10 = sum((1-Y)*Z)/sum(1-Y), z1y11 = sum(Y*Z)/sum(Y)) -> samplemean
  MSE_sample = apply((truemean[,-1] - samplemean[,-1])^2, MARGIN = 2, FUN = mean, na.rm =T)
  
  
  # univariate model (uncomment if necessary)
  if(useuni == T)
  {
    um_z = lmer(as.formula(fz_uni), data = sampledata)
    um_y1 = glmer(Y ~ X2 + (1|group), data = sampledata, family = binomial, nAGQ = 20)

    zpredict = predict(um_z, unsampledata)
    By_hat_uni <- fixef(um_y1) %>% unlist()

    Vy_uni <- rep(0, 40)
    Vy_uni[rownames(ranef(um_y1)$group ) %>% as.numeric()] <- ranef(um_y1)$group %>% unlist
    SigmaV_hat <- var(Vy_uni)

    EBPfun_uni(sampledata, unsampledata, By_hat_uni, SigmaV_hat, m, 50, zpredict) %>%
     apply(., MARGIN = c(1,2), FUN = mean, na.rm = T) %>% as.vector -> unimodelmean
    MSE_uni = apply((truemean[,-1] - unimodelmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  }

  # Joint Model (our approach): 
  rcppEBPfun(sampledata, unsampledata, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, m = 200, t)  -> temp
  jointmean = apply(temp, MARGIN = c(1,2), FUN = mean, na.rm = T)  %>% as.vector
  jointvar = apply(temp, MARGIN = c(1,2), FUN = var, na.rm = T)  %>% as.vector

  MSE_joint = apply((truemean[,-1] - jointmean)^2, MARGIN = 2, FUN = mean, na.rm =T)

  # compute MSE with bootstrap method if needed: 
  # (NOTE: use for loop in runME_oneZoneY if you need to use Bootstrap method so that you can do parallel here to save computation time)
  
  # for computing MSE and BR 
  if(useuni == T) 
  {
    # result = cbind(MSE_sample, MSE_uni, MSE_joint) 
    temp = temp %>% unlist %>% matrix(., nrow = D * 5)
    result = cbind(truemean[,-1] %>% unlist, samplemean[,-1] %>% unlist, unimodelmean, jointmean, jointvar)            
  }       
  if(useuni == F)
  {
    # result = cbind(MSE_sample, MSE_joint) 
    temp = temp %>% unlist %>% matrix(., nrow = D * 5)
    result = cbind(truemean[,-1] %>% unlist, samplemean[,-1] %>% unlist, jointmean, jointvar) 
  }

  if(bootstrap == T)
  {
    mclapply(1:b, FUN = function(x) oneboot(sampledata, unsampledata, Bz_hat, By_hat, 
                                         SigmaZ_hat, P_hat, SigmaUV_hat, m, t, EMsize), 
             mc.cores = cores) %>% do.call(cbind,.)-> boot

    saveRDS(cbind(truemean[,-1] %>% unlist, jointmean, jointvar, boot), 
            file = paste("result/bootstrapMSE_", x, ".rds", sep = ""))
  }

}
  
if(doEBP == FALSE)
  return(rbind(initialest, c(Bz_hat, By_hat, P_hat, SigmaZ_hat, as.vector(SigmaUV_hat))))
if(doEBP == TRUE)
  return(result)
}

