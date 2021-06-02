# compute the log density for sample data in each group
# including f_yz, f_z, and f_y 
groupdenscpp <-'
int m = Rcpp::as<int>(m0);            // MC replicates 
int D = Rcpp::as<int>(D0);            // group size 
double P = Rcpp::as<double>(P0);
double sigmaZ = Rcpp::as<double>(sigmaZ0);            
arma::vec Z = Rcpp::as<arma::vec>(Z0); 
arma::vec Y = Rcpp::as<arma::vec>(Y0); 
arma::cube ub = Rcpp::as<arma::cube>(rand_matrix);    // generate random effect  D * m * 2
arma::vec delta = Rcpp::as<arma::vec>(delta_fix0);    // x * beta_z
arma::vec mu = Rcpp::as<arma::vec>(mu_fix0);          // x * beta_y
arma::vec id = Rcpp::as<arma::vec>(id0);              // group id 
arma::vec subid = Rcpp::as<arma::vec>(subid0);        // subgroup id (1:5: H1,K1,H0,K0, M)

int n = Z.n_rows;                  // sample size
double gx; 
double u; 
double b;
double logdens;
arma::mat result = arma::zeros<arma::mat>(m, D);     //result table ???

for(int i = 1; i<=m; i++)         
{
  for(int j = 1; j<=n; j++)
  {
  u = ub(id(j-1)-1, i-1, 0);      // random effect for Z
  b = ub(id(j-1)-1, i-1, 1);      // random effect for Y
  gx =  P * (delta(j-1) + u) + mu(j-1) + b + sigmaZ * pow(P, 2)/2;    // g(x) 
  logdens = Y(j-1) * gx -log(1+exp(gx));    
  if(subid(j-1) < 3)  // H, K
  {
  logdens += -0.5 * log(2*PI*sigmaZ) - pow(Z(j-1) - delta(j-1) - u - Y(j-1)*sigmaZ*P, 2)/sigmaZ/2; 
  }
  if(subid(j-1) == 2)  // K1
  {
  logdens += log(1+exp(Z(j-1)*P + mu(j-1) + b));
  }
  result(i-1, id(j-1)-1) += logdens;
  } 
}
return Rcpp::wrap(result);
'

rcppgroupdens <- cxxfunction(signature(m0 = "int", D0 = "int", P0 = "numeric", sigmaZ0 = "numeric", 
                                       Z0 = "numeric", Y0 = "numeric", rand_matrix = "numeric",
                                       delta_fix0 = "numeric", mu_fix0 = "numeric", id0 = "int", subid0 = "int"),
                             groupdenscpp,plugin="RcppArmadillo", verbose = TRUE)

SigmaUVweightavgcpp <- ' 
arma::mat xrand = Rcpp::as<arma::mat>(rand_origin);  // mD * 2
arma::vec xweight = Rcpp::as<arma::vec>(weight);    // D*m
int n = xweight.n_rows;
arma::mat result = arma::zeros<arma::mat>(2,2);
for(int i = 1; i<=2; i++)
for(int j = 1; j<=2; j++)
for(int k = 1; k<=n; k++)
result(i-1,j-1) += xrand(k-1,i-1) * xrand(k-1,j-1) * xweight(k-1);
return Rcpp::wrap(result);
'

rcppSigmaUVweightavg <- cxxfunction(signature(rand_origin="numeric", weight = "numeric"),
                                    SigmaUVweightavgcpp,plugin="RcppArmadillo", verbose = TRUE)

# Ydraw n * 3 (first column denote the label 0 or 1, second column denote the Yid 1 or 2, third column is condition mean of Z)
YZdrawcpp <- ' 
  arma::vec xdelta = Rcpp::as<arma::vec>(delta_new);
  arma::vec xmu = Rcpp::as<arma::vec>(mu_new);
  arma::mat xprob = Rcpp::as<arma::mat>(Yprob_new);
  double xP = Rcpp::as<double>(P_hat);
  double xsimgaZ = Rcpp::as<double>(SigmaZ_hat);
  int n = xmu.n_rows;
  arma::mat Ydraw = arma::zeros<arma::mat>(n,3);
  arma::vec rand = Rcpp::as<arma::vec>(rand_vec);
  for(int i = 1; i<=n; i++)
  {
    if(rand(i-1) < xprob(i-1, 1))
    {
      Ydraw(i-1,0) = 1;
      Ydraw(i-1,1) = 2;
    }
    else
      Ydraw(i-1,1) = 1;
    Ydraw(i-1,2) = xdelta(i-1) + Ydraw(i-1,0) * xP * xsimgaZ;
  }
  return Rcpp::wrap(Ydraw);
'

rcppYZdraw <- cxxfunction(signature(delta_new = "numeric", mu_new = "numeric",
                                    Yprob_new = "numeric", P_hat = "numeric", SigmaZ_hat = "numeric", rand_vec = "numeric"), 
                          YZdrawcpp, plugin = "RcppArmadillo", verbose = TRUE)

groupdens_expectation_cpp <-'
double P = Rcpp::as<double>(P0);
double sigmaZ = Rcpp::as<double>(sigmaZ0);   
arma::vec Z = Rcpp::as<arma::vec>(Z0); 
arma::vec Y = Rcpp::as<arma::vec>(Y0); 
arma::cube ub = Rcpp::as<arma::cube>(rand_matrix_posterior);    //  m * D * 2
arma::vec delta = Rcpp::as<arma::vec>(delta_fix0);    // x * beta_z
arma::vec mu = Rcpp::as<arma::vec>(mu_fix0);          // x * beta_y
arma::vec id = Rcpp::as<arma::vec>(id0);              // group id 
arma::vec subid = Rcpp::as<arma::vec>(subid0);        // subgroup id (1:5: H1,K1,H0,K0, M)

int n = Z.n_rows;                  // sample size
int m = ub.n_rows;            // MC replicates 
double gx; 
double u; 
double b;
double logdens;               //logdens for one sample                  
double result;            // final output

result = 0; 
for(int j = 1; j<=n; j++)        // each sample 
{
  logdens = 0;                   // first define logdens (later random part)
  for(int i = 1; i<=m; i++)         //each MC draw 
  {
  u = ub(i-1, id(j-1)-1, 0);      // random effect for Z
  b = ub(i-1, id(j-1)-1, 1);      // random effect for Y
  gx =  P * (delta(j-1) + u) + mu(j-1) + b + sigmaZ * pow(P, 2)/2;    // g(x)  
  logdens += -log(1+exp(gx)) + Y(j-1) * gx;     // include first term random part
  if(subid(j-1) < 3) 
  {
    logdens += -0.5 * log(2*PI*sigmaZ) - pow(Z(j-1) - delta(j-1) - Y(j-1)*P*sigmaZ - u, 2)/sigmaZ/2;   // random part
  }
  if(subid(j-1) == 2)  // K
  {
    logdens += log(1 + exp(Z(j-1) * P + mu(j-1) + b));  
  }
  } 
  result -= logdens/m;   //expectation of negative loglikehood
}

return Rcpp::wrap(result);
'

rcppgroupdens_expectation_cpp <- cxxfunction(signature(P0 = "numeric", sigmaZ0 = "numeric", Z0 = "numeric", Y0 = "numeric", 
                                                       rand_matrix_posterior = "numeric", delta_fix0 = "numeric", 
                                                       mu_fix0 = "numeric", id0 = 'int', subid0 = "int"),
                             groupdens_expectation_cpp,plugin="RcppArmadillo", verbose = TRUE)




groupdens_dev1_cpp <-'
double P = Rcpp::as<double>(P0);
double sigmaZ = Rcpp::as<double>(sigmaZ0);   
int num_Zcoef = Rcpp::as<int>(num_Zcoef0);  
arma::vec Z = Rcpp::as<arma::vec>(Z0); 
arma::vec Y = Rcpp::as<arma::vec>(Y0); 
arma::mat X = Rcpp::as<arma::mat>(X0);      //  matrix   
arma::cube ub = Rcpp::as<arma::cube>(rand_matrix_posterior);    // m * D* 2
arma::vec delta = Rcpp::as<arma::vec>(delta_fix0);    // x * beta_z
arma::vec mu = Rcpp::as<arma::vec>(mu_fix0);          // x * beta_y
arma::vec id = Rcpp::as<arma::vec>(id0);              // group id 
arma::vec subid = Rcpp::as<arma::vec>(subid0);        // subgroup id (1:3: H, K, M)

int n = Z.n_rows;                  // sample size
int m = ub.n_rows;            // MC replicates 
double gx; 
double u; 
double b;
double dbetaz0 = 0;
arma::vec dbetaz1  = arma::zeros<arma::vec>(num_Zcoef);      // vector 

double dbetay0 = 0; 
double dbetay1 = 0; 
double dP = 0; 
double dsigmaZ = 0; 
double partA;
double partB; 

arma::vec result = arma::zeros<arma::vec>(4+num_Zcoef+1);      // betaz, betay, P, sigmaZ;

for(int j = 1; j<=n; j++)        // each sample 
{
  for(int i = 1; i<=m; i++)         //each MC draw 
  {
    u = ub(i-1, id(j-1)-1, 0);      // random effect for Z
    b = ub(i-1, id(j-1)-1, 1);      // random effect for Y
    gx =  P * (delta(j-1) + u) + mu(j-1) + b + sigmaZ * pow(P, 2)/2;    // g(x) 
    partB = -exp(gx)/(1+exp(gx));    
    dbetaz0 += (Y(j-1) + partB) * P;
       for(int k = 1; k <= num_Zcoef; k++)
       {
            dbetaz1(k-1) += (Y(j-1) + partB) * P * X(j-1, k-1); 
       }
    dbetay0 += Y(j-1) + partB; 
    dbetay1 += (Y(j-1) + partB) * X(j-1, 0); 
    dP += (Y(j-1) + partB) * (delta(j-1) + u + sigmaZ* P);   // add random part here 
    dsigmaZ += (Y(j-1) + partB) * pow(P, 2)/2; 
  if(subid(j-1) < 3)
  { 
    partA = (Z(j-1) - delta(j-1) - u - Y(j-1) * P *sigmaZ)/sigmaZ; 
    dbetaz0 += partA;
    for(int k = 1; k <= num_Zcoef; k++)
    {
      dbetaz1(k-1) += partA * X(j-1, k-1);
    }
    dP += partA * Y(j-1) * sigmaZ; 
    dsigmaZ += -0.5/sigmaZ + pow(partA, 2)/2 + partA * Y(j-1)*P; 
  }
  if(subid(j-1) == 2)  // K
  {
    partB = exp(Z(j-1) * P + mu(j-1) + b)/(1 + exp(Z(j-1) * P + mu(j-1) + b));
    dbetay0 += partB;
    dbetay1 += partB * X(j-1,0); 
    dP += partB * Z(j-1);
  }
  } 
  result(0) -= dbetaz0/m;
  for(int k = 1; k<= num_Zcoef; k++)
  {
    result(k) -= dbetaz1(k-1)/m;
  }
  result(num_Zcoef+1) -= dbetay0/m;
  result(num_Zcoef+2) -= dbetay1/m; 
  result(num_Zcoef+3) -= dP/m;
  result(num_Zcoef+4) -= dsigmaZ/m;   //derivative of negative loglikelihood
}

return Rcpp::wrap(result);
'

rcppgroupdens_dev1_cpp  <- cxxfunction(signature(P0 = "numeric", sigmaZ0 = "numeric", num_Zcoef0 = 'int',
                                                 Z0 = "numeric", Y0 = "numeric", X0 = "numeric", rand_matrix_posterior = "numeric",
                                                 delta_fix0 = "numeric", mu_fix0 = "numeric", id0 = "int", subid0 = "int"),
                                       groupdens_dev1_cpp ,plugin="RcppArmadillo", verbose = TRUE)

rcppfulldevfun <- function(data, initparam, rand_matrix_posterior)
{
  bZinit = initparam[1:(1+num_Zcoef)]      
  bYinit = initparam[(2+num_Zcoef):(2+num_Zcoef+num_Ycoef)]
  Pinit = initparam[3+num_Zcoef+num_Ycoef]
  SigmaZinit = initparam[4+num_Zcoef+num_Ycoef]
  delta_fix = as.matrix(cbind(1, data[,2:(1+num_Zcoef)])) %*% bZinit
  mu_fix = cbind(1, data$X2) %*% bYinit
    X0 = data[, 2:(1+num_Zcoef)] %>% as.matrix()
  dev = rcppgroupdens_dev1_cpp(P0 = P, sigmaZ0 = SigmaZ, Z0 = data$Z, Y0 = data$Y, 
                               X0 = X0, rand_matrix_posterior = rand_matrix_posterior, delta_fix0 = delta_fix, 
                               mu_fix0 = mu_fix, id0 = data$group, subid0 = data$subgroupid)
  return(dev)
}

rcpploglikefun <- function(data, initparam, rand_matrix_posterior)   #(check)
{
  bZinit = initparam[1:(1+num_Zcoef)]    
  bYinit = initparam[(2+num_Zcoef):(2+num_Zcoef+num_Ycoef)]
  Pinit = initparam[3+num_Zcoef+num_Ycoef]
  SigmaZinit = initparam[4+num_Zcoef+num_Ycoef]
  delta_fix = as.matrix(cbind(1, data[,2:(1+num_Zcoef)])) %*% bZinit
  mu_fix = cbind(1, data$X2) %*% bYinit
  r1 = rcppgroupdens_expectation_cpp(P0 = Pinit, sigmaZ0 = SigmaZinit, Z0 = data$Z, Y0 = data$Y, 
                                     rand_matrix_posterior = rand_matrix_posterior, delta_fix0 = delta_fix, mu_fix0 = mu_fix, 
                                     id0 = data$group, subid0 = data$subgroupid)
  return(r1)
}

rcppupdatefix <- function(data, Bz_hat, By_hat, SigmaZ_hat, P_hat, SigmaUV_hat, maxiter, m)
{
  rand_matrix <- rmvnorm(m * D, mean = rep(0, 2), sigma = SigmaUV_hat) %>% array(., dim = c(D, m, 2)) # 40 * 10000 
  
  delta_fix = as.matrix(cbind(1, data[, 2:(1+num_Zcoef)])) %*% Bz_hat
  mu_fix = cbind(1, data$X2) %*% By_hat

  # compute the posterior weight 
  logdens <- rcppgroupdens(m0 = m, D0 = D, P0 = P_hat, sigmaZ0 = SigmaZ_hat, Z0 = data$Z, Y0 = data$Y, 
                           rand_matrix = rand_matrix, delta_fix0 = delta_fix, mu_fix0 = mu_fix, 
                           id0 = data$group, subid0 = data$subgroupid)
  quant_logdens = apply(logdens, MARGIN = 2, FUN = function(x) quantile(x, 0.7)) 
  dens = sweep(logdens, MARGIN = 2, FUN = "-", STATS = quant_logdens) %>% exp()
  dens[is.infinite(dens)] <- 0
  weight = sweep(dens, STATS= colSums(dens), FUN = "/" , MARGIN = 2)    # posterior weight
  # next generate posterior samples from posterior distribution
  posterior_ID <- sapply(1:D, FUN = function(x) base::sample(1:dim(weight)[1], dim(weight)[1], replace = TRUE, prob = weight[,x])) %>% 
    as.vector() %>% cbind(rep(1:D, each = dim(weight)[1]), .)  
  rand_matrix_posterior <- c(rand_matrix[,,1][posterior_ID], rand_matrix[,,2][posterior_ID]) %>% array(., dim = c(m, D, 2)) # 1000 * 40 * 2
  
  initparam <- c(Bz_hat, By_hat, P_hat, SigmaZ_hat)
  
  if(usedev == TRUE)
    o1 <- optim(par = initparam, fn = rcpploglikefun, data = data, method = "BFGS", gr = rcppfulldevfun,
                rand_matrix_posterior = rand_matrix_posterior, control = list(maxit = maxiter))
  if(usedev == FALSE)
    o1 <- optim(par = initparam, fn = rcpploglikefun, data = data, 
                rand_matrix_posterior = rand_matrix_posterior, control = list(maxit = maxiter))

  return(o1$par)  #(bz, by, P, sigmaZ)
}






SigmaUVweightavgcpp <- ' 
arma::mat xrand = Rcpp::as<arma::mat>(rand_origin);
arma::vec xweight = Rcpp::as<arma::vec>(weight);
int n = xweight.n_rows;
arma::mat result = arma::zeros<arma::mat>(2,2);
for(int i = 1; i<=2; i++)
for(int j = 1; j<=2; j++)
for(int k = 1; k<=n; k++)
result(i-1,j-1) += xrand(k-1,i-1) * xrand(k-1,j-1) * xweight(k-1);
return Rcpp::wrap(result);
'

rcppSigmaUVweightavg <- cxxfunction(signature(rand_origin="numeric", weight = "numeric"),
                                    SigmaUVweightavgcpp,plugin="RcppArmadillo", verbose = TRUE)

# this function uses two rcpp functions:
rcppSigmaUVfun <- function(data, SigmaZ_hat, P_hat, Bz_hat, By_hat, SigmaUV_hat, m)
{
  rand_origin = rmvnorm(n = m*D, mean = rep(0, 2), sigma = SigmaUV_hat)   
  rand_matrix  = rand_origin %>% array(., dim = c(D, m, 2))  # D*m*2

  delta_fix = as.matrix(cbind(1, data[, 2:(1+num_Zcoef)])) %*% Bz_hat  
  mu_fix = cbind(1, data$X2) %*% By_hat    # now is a vector

  rcppgroupdens(m0 = m, D0 = D, P0 = P, sigmaZ0 = SigmaZ_hat, Z0 = data$Z, Y0 = data$Y, 
                rand_matrix = rand_matrix, delta_fix0 = delta_fix, 
                mu_fix0 = mu_fix, id0 = data$group, subid0 = data$subgroupid) -> logdens  # m * D
  rm(rand_matrix)  
  quant_logdens = apply(logdens, MARGIN = 2, FUN = function(x) quantile(x, 0.7)) 
  dens = sweep(logdens, MARGIN = 2, FUN = "-", STATS = quant_logdens) %>% exp()
  dens[is.infinite(dens)] <- 0
  weight = t(sweep(dens, STATS= colSums(dens), FUN = "/" , MARGIN = 2)/D)   # D * m 
  
  return(rcppSigmaUVweightavg(rand_origin, as.vector(weight)))
}


# Ydraw n * 3 (first column denote the label 0 or 1, second column denote the Yid 1 or 2, third column is condition mean of Z)
YZdrawcpp <- ' 
  arma::vec xdelta = Rcpp::as<arma::vec>(delta_new);
  arma::vec xmu = Rcpp::as<arma::vec>(mu_new);
  arma::mat xprob = Rcpp::as<arma::mat>(Yprob_new);
  double xP = Rcpp::as<double>(P_hat);
  double xsimgaZ = Rcpp::as<double>(SigmaZ_hat);
  int n = xmu.n_rows;
  arma::mat Ydraw = arma::zeros<arma::mat>(n,3);
  arma::vec rand = Rcpp::as<arma::vec>(rand_vec);
  for(int i = 1; i<=n; i++)
  {
    if(rand(i-1) < xprob(i-1, 1))
    {
      Ydraw(i-1,0) = 1;
      Ydraw(i-1,1) = 2;
    }
    else
      Ydraw(i-1,1) = 1;
    Ydraw(i-1,2) = xdelta(i-1) + Ydraw(i-1,0) * xP * xsimgaZ;
  }
  return Rcpp::wrap(Ydraw);
'

rcppYZdraw <- cxxfunction(signature(delta_new = "numeric", mu_new = "numeric",
                                    Yprob_new = "numeric", P_hat = "numeric", SigmaZ_hat = "numeric", rand_vec = "numeric"), 
                          YZdrawcpp, plugin = "RcppArmadillo", verbose = TRUE)