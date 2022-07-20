## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  n <- 40
#  simp <- rmultinom(n,1,c(0.6,0.4))
#  lab <- as.factor(apply(t(simp),1,which.max))
#  
#  
#  #parameter comes from multinomial
#  p <- 11
#  mu1 <- c(-2.8,-1.3,-1.6,-3.9,-2.6,-2.9,-2.5,-2.7,-3.1,-2.9)
#  B1 <- matrix(c(1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0),nrow = p-1)
#  T1 <- diag(c(2.9,0.5,1))
#  D1 <- diag(c(0.52, 1.53, 0.56, 0.19, 1.32, 1.77, 0.6, 0.53, 0.37, 0.4))
#  cov1 <- B1%*%T1%*%t(B1)+D1
#  
#  mu2 <- c(1.5,-2.7,-1.1,-0.4,-1.4,-2.6,-3,-3.9,-2.7,-3)
#  B2 <- matrix(c(1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,0,0,0,0),nrow = p-1)
#  T2 <- diag(c(0.2,0.003,0.15))
#  D2 <- diag(c(0.01, 0.62, 0.45, 0.01, 0.37, 0.42, 0.08, 0.16, 0.23, 0.27))
#  cov2 <- B2%*%T2%*%t(B2)+D2
#  
#  df <- matrix(0,nrow=n,ncol=p-1)
#  for (i in 1:n) {
#    if(lab[i]==1){df[i,] <- rmvnorm(1,mu1,sigma = cov1)}
#    else if(lab[i]==2){df[i,] <- rmvnorm(1,mu2,sigma = cov2)}
#  }
#  
#  
#  
#  f_df <- cbind(df,0)
#  z <- exp(f_df)/rowSums(exp(f_df))
#  
#  
#  W_count <- matrix(0,nrow=n,ncol=p)
#  for (i in 1:n) {
#    W_count[i,] <- rmultinom(1,runif(1,10000,20000),z[i,])
#  
#  }
#  

## ----eval=FALSE---------------------------------------------------------------
#  range_G <- 2 #define the number of components
#  range_Q <- 2 #define the possible number of bicluster.
#  cov_str <- "GUU" #select the model you want to fit
#  #It will fit GUU model with G=2, Q=c(2,2)
#  res <- lnmbiclust(W_count=W_count, range_G=range_G, range_Q=range_Q, model=cov_str)
#  #where res will be a list contain all parameters.
#  

## ----eval=FALSE---------------------------------------------------------------
#  range_G <- c(2:3)
#  range_Q <- c(2:3)
#  cov_str <- c("UUU", "GGC")
#  res <- lnmbiclust(W_count=W_count, range_G=range_G, range_Q=range_Q, model=cov_str, criteria="BIC")
#  best_model=res$best_model
#  

## ----eval=FALSE---------------------------------------------------------------
#  range_G <- 2
#  range_Q <- c(2:3)
#  cov_str <- "UUU"
#  res <- lnmbiclust(W_count=W_count, range_G=range_G, range_Q=range_Q, model=cov_str, criteria="BIC",permutation=TRUE)
#  res$best_model
#  

## ----eval=FALSE---------------------------------------------------------------
#  res$all_fitted_model
#  

## ----eval=FALSE---------------------------------------------------------------
#  range_G <- c(2:3)
#  range_Q <- c(2:3)
#  cov_str <- c("UUU", "GUC")
#  res <- lnmfa(W_count=W_count, range_G=range_G, range_Q=range_Q, model=cov_str, criteria="BIC")
#  best_model=res$best_model
#  model_output=res$all_fitted_model
#  

## ----eval=FALSE---------------------------------------------------------------
#  range_G <- c(2:3)
#  range_tunning=seq(0.5,0.7,length.out=10)
#  range_Q <- 2
#  cov_str <- "UUU"
#  res <- plnmfa(W_count=W_count, range_G=range_G, range_Q=range_Q, model=cov_str, criteria="BIC",
#                range_tuning = range_tuning)
#  best_model=res$best_model
#  model_output=res$all_fitted_model
#  

