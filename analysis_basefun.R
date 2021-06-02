AvMSE_BR <- function(result, pos = 3)
{
  S = length(result)   # list
  MSE_MC = rep(0, D * 5)
  Bias_MC = rep(0, D * 5)
  NArecord = rep(0, D * 5)
  for(i in 1:S)
  {
    temp = result[[i]]
    NArecord = NArecord + as.numeric((is.na(temp[,1]) + is.na(temp[,pos])) > 0)
    temp[is.na(temp)] <- 0   # fill in 0 
    
    MSE_MC = MSE_MC + (temp[,1] - temp[,pos])^2 
    Bias_MC = Bias_MC + (temp[,1] - temp[,pos])
  }
  MSE_MC = matrix(MSE_MC/(rep(S, D*5) - NArecord), nrow = D)
  Bias_MC = matrix(Bias_MC/(rep(S, D*5) - NArecord), nrow = D)
  AvMSE = colMeans(MSE_MC, na.rm = T)
  BR = colMeans(Bias_MC^2, na.rm = T)/AvMSE
  return(cbind(AvMSE * 100, BR * 100) %>% "rownames<-"(c("Z", "Y0", "Y1", "Z|Y0", "Z|Y1")) %>% 
      "colnames<-"(c("MSE*100", "BR*100%")))
}


relative_bias <- function(S, D_each)  # D_each = list(1:5, 6:10, 11:15, 16:20)
{
  MSE_MC = rep(0, D * 5)
  MSE_est = rep(0, D * 5)
  NArecord = rep(0, D * 5)
  for(i in 1:S)
  {
    temp = readRDS(paste("result/bootstrapMSE_", i, ".rds", sep = ""))  # dim = (D, 3 + 2*b)
    NArecord = NArecord + as.numeric((is.na(temp[,1]) + is.na(temp[,2]) + is.na(temp[,3])) > 0)
    temp1 = temp[,1:3]
    temp1[is.na(temp1)] <- 0
    temp[,1:3] <- temp1
    MSE_MC = MSE_MC + (temp[,1] - temp[,2])^2 
    
    Varpart = temp[,3]
    bias2 = ((temp[,seq(4, length = b, by = 3)] - temp[,2])^2 ) %>% rowMeans(na.rm = T)
    # Varbias = temp[,seq(5, length = b, by = 3)] %>% rowMeans(na.rm = T) %>% "-"(Varpart)
    MSE_est = MSE_est + Varpart + bias2  
  }
  MSE_MC = matrix(MSE_MC/(rep(S, D*5) - NArecord), nrow = D)
  MSE_est = matrix(MSE_est/rep(S, D*5), nrow = D)
  
  RB = NULL
  for(i in 1:length(D_each))
  {
    MSE_MC_g = MSE_MC[D_each[[i]], ]
    MSE_est_g = MSE_est[D_each[[i]], ]
    RB_group = colMeans(MSE_est_g - MSE_MC_g, na.rm = T)/colMeans(MSE_MC_g) * 100
    RB = rbind(RB, RB_group)
  }
  
  RB %>% "rownames<-"(c("100", "200", "300", "400")) %>% return()
}
