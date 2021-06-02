# Stan model 
model_SAE = "
data {
  int<lower=0> D;                 // num of groups  
  int<lower=0> n;                 // sample size
  int<lower=0> groupid[n];        // group index the sample unit belongs to
  int<lower=0, upper=1> Y[n];     // observed Y
  vector[n] Z1;                   // observed Z
  vector[n] X;                    // auxillary variable 
  matrix[2,2] Sigma0;             // prior for inverse wishart 
  real<lower=0> b0;               // 
}
parameters {
  real beta_z1[2];                // fixed coefficient parameters 
  real beta_y[2]; 
  vector[2] U[D];                 // random effect U = (u1, u2), u1 for Z1, u2 for Z2
  cov_matrix[2] preU;             // precision matrix for U; 
  real<lower=0> varE1;            // variance for Z1
  real<lower=-1, upper=1> rhoE;   // correlation between E1 and E2 
}
transformed parameters {
  vector[n] u1all;                // u1 vector
  vector[n] u2all; 
  vector[n] mu2;               // conditional expectation for latent variable Z2
  for (i in 1:n) {
    u1all[i] = U[groupid[i]][1]; 
    u2all[i] = U[groupid[i]][2]; 
  }
  mu2 = beta_y[1] + beta_y[2] * X + u2all + rhoE/sqrt(varE1) * (Z1 - beta_z1[1] - beta_z1[2] * X - u1all);   

}
model {
  preU ~ inv_wishart(b0, b0 * Sigma0);    //precision has an inverse wishart distribution
  varE1 ~ inv_gamma(1, 1);                //variance E1 has a inverse gamma prior
  U ~ multi_normal_prec([0, 0], preU);    //generate random effect U
  Z1 ~ normal(beta_z1[1] + beta_z1[2] * X + u1all, sqrt(varE1));  // generate Z1
  Y  ~ bernoulli(Phi(mu2/sqrt(1-rhoE^2)));         // generate Y 
}
"

stan_model1 = stan_model(model_code = model_SAE)


# Cpp code: 
BayesrcppEBPonedraw <-function(data, fitting, Bz_hat, By_hat, SigmaZ_hat, rhoE_hat, deltanew_fix, munew_fix)
{
  n = length(deltanew_fix)
  VarE = matrix(c(SigmaZ_hat, rep(rhoE_hat * sqrt(SigmaZ_hat), 2), 1), nrow = 2)
  
  E = rmvnorm(n, mean = c(0, 0), sigma = VarE)
  Z_draw = deltanew_fix + E[,1]
  Y_draw = as.numeric((munew_fix + E[,2]) > 0)
  data.frame(fitting, Y = Y_draw, Z = Z_draw, Yid = Y_draw %>% as.factor) %>% 
    rbind(., data.frame(data[,1:5])) %>% 
    group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y), y11mean = mean(Y), 
                                  z1y10 = sum((1-Y)*Z)/sum(1-Y), z1y11 = sum(Y*Z)/sum(Y)) %>% as.matrix() -> mean_draw
  return(mean_draw[,-1])
}


BayesrcppEBPfun <- function(data, fitting, bayesest)
{
  Bz_hat <- bayesest[1:2]
  By_hat <- bayesest[3:4]
  SigmaZ_hat <- bayesest[5]
  rhoE_hat <- bayesest[6]
  randef <- bayesest[7:(6+2*D)] %>% matrix(., nrow = D)   #u1 for z, u2 for y   (D * 2) 
  
  deltanew1_fix = cbind(1,fitting$X2) %*% Bz_hat + randef[,1][fitting$group]
  deltanew2_fix = cbind(1,fitting$X2) %*% By_hat + randef[,2][fitting$group]
  
  temp = BayesrcppEBPonedraw(data, fitting, Bz_hat, By_hat, SigmaZ_hat, rhoE_hat, deltanew1_fix, deltanew2_fix)
  return(temp)
}


Bayes_onesimu <- function(x)
{
  if(noninfor == FALSE)
  {
    sampledata = readRDS(paste("data/SAE_sample_infor_", x, ".rds", sep = ""))
  }
  if(noninfor == TRUE)
  {
    sampledata = readRDS(paste("data/SAE_sample_", x, ".rds", sep = ""))
  }
  
  popdata = readRDS(paste("data/SAE_pop_", x, ".rds", sep = ""))
  unsampledata = anti_join(popdata, sampledata)[,1:2]

  stan_data1 = list(D = D, n = dim(sampledata)[1], groupid = sampledata$group, 
                    Y = sampledata$Y, Z1 = sampledata$Z, X = sampledata$X2, Sigma0 = diag(c(1, 1)), 
                    b0 = 3)
  
  stan_r = sampling(stan_model1, stan_data1, c("beta_z1", "beta_y", "varE1", "rhoE", "U"), chains = 4, iter = 10000)
  samples = stan_r@sim$samples
  posterior <- NULL
  for(j in 1:4)
  {
    posterior <- rbind(posterior, sapply(1:206, FUN = function(k) samples[[j]][[k]][5001:10000]))
  }
  saveRDS(posterior, file = paste("rstan/posterior_", x, ".rds", sep = ""))  # save posterior result 
  
  # true group mean 
  popdata %>% group_by(group) %>% summarise(z1mean = mean(Z), y10mean = mean(1-Y), y11mean = mean(Y), 
                                            z1y10 = sum((1-Y)*Z)/sum(1-Y), z1y11 = sum(Y*Z)/sum(Y)) -> truemean
  
  mclapply(1:20000, FUN = function(x) BayesrcppEBPfun(sampledata, unsampledata, posterior[x,]),
           mc.cores = cores) %>% unlist() %>% array(., dim = c(D, 5, 20000)) -> temp 
  
  bayesmean = apply(temp, MARGIN = c(1,2), FUN = mean) %>% unlist() %>% as.vector()
  
  MSE_bayes = apply((truemean[,-1] - bayesmean)^2, MARGIN = 2, FUN = mean, na.rm =T)
  return(cbind(truemean[,-1] %>% unlist, bayesmean))
}




