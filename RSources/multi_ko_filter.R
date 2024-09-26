
### W rank for multi knockoff
# Xk.list: list
ko.filter.multiko <- function(X, Xk.list, y, q, method = "lasso.lcd", gap = median, sparsity = 0.1, data.log = T){
  n = nrow(X)
  p = ncol(X)
  
  
  
  if("list" %in% class(Xk.list)){
    Xk = do.call("cbind", Xk.list) 
    M.ko = length(Xk.list)
  } 
  if("matrix" %in% class(Xk.list) | "dataframe" %in% class(Xk.list)){
    Xk = Xk.list 
    M.ko = ncol(Xk)/ncol(X)
  } 

  
  if(method == "lasso.lcd"){
    coefs = cv_coeffs_glmnet(cbind(X, Xk), y, family='gaussian', parallel = F)
    Z <- coefs[2:((1+M.ko)*p+1)]
  }
  else if(method == "lasso.logistic"){
    coefs = cv_coeffs_glmnet(cbind(X, Xk), y, family='binomial', parallel = F)
    Z <- coefs[2:((1+M.ko)*p+1)]
  }
  else if(method == "p"){
    raw.p <- apply(cbind(X, Xk),2, function(Xj){
      p <- summary(lm(Xj~y))[["coefficients"]][2,4]
    })
    Z <- -log10(raw.p)
  }
  
  Z.matrix <- matrix(data = Z, nrow = 1+M.ko, ncol = p, byrow = T)
  Z.matrix.abs <- abs(Z.matrix)
  
  # kappa denote the index of the original (denoted as 0) or knockoff feature that has the largest importance score
  kap <- apply(Z.matrix.abs, 2, which.max)-1
  
  # if orginal var has the largest importance value, I.max = 1, 0 otherwise
  I.max <- ifelse(kap==0, yes = 1, no = 0)
  # Median of knockoffs' importance value
  T_m_median <- apply(Z.matrix.abs[2:(M.ko+1),], 2, median)
  
  # tau denote the difference between the largest importance score and the median of the remaining importance scores.
  tau <- rep(0, p)
  for(i in 1:p){
    tau[i] <- Z.matrix.abs[kap[i]+1, i] - gap(Z.matrix.abs[-(kap[i]+1), i])
  }

  W = tau * I.max
  
  # find the knockoff threshold T
  t = sort(c(0, tau))

  ratio = c(rep(0, p))

  for (j in 1:p) {
    ratio[j] = (1/M.ko+1/M.ko*sum((tau >= t[j]) & (kap>=1)))/max(1, sum((tau >= t[j])& (kap == 0)))
    #ratio[j] = (1/M.ko*sum((tau >= t[j]) & (kap>=1)))/max(1, sum((tau >= t[j])& (kap == 0)))
  }
  
  id = which(ratio <= q)[1]

  if(length(id) == 0){
    T = Inf
  } else {
    T = t[id]
  }
  
  # set of discovered variables
  S = integer0_test(which(W > T))

  return(list(S = S, W = W, tau = tau, kap = kap, T =T, Z.matrix = Z.matrix))
}

# 
# cv_coeffs_glmnet_mod <- function(X, y, nlambda=500, intercept=T, parallel=T, family = "gaussian") {
#   # Standardize variables
#   X = scale(X)
#   
#   n = nrow(X); p = ncol(X)
#   
#   if (!methods::hasArg(lambda) ) {
#     if( identical(family, "gaussian") ) {
#       if(!is.numeric(y)) {
#         stop('Input y must be numeric.')
#       }
#       # Unless a lambda sequence is provided by the user, generate it
#       lambda_max = max(abs(t(X) %*% y)) / n
#       lambda_min = lambda_max / 2e3
#       k = (0:(nlambda-1)) / nlambda
#       lambda = lambda_max * (lambda_min/lambda_max)^k
#     }
#     else {
#       lambda = NULL
#     }
#   }
#   
#   cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, intercept=intercept,
#                                      standardize=F,standardize.response=F, parallel=parallel)
#   
#   coef(cv.glmnet.fit, s = "lambda.min")
# }

# gogohugo cheapknockoff filtering ======
library(CVXR)

source("../externalCodes/cheapknockoff-master/R/generate_X.R")
source("../externalCodes/cheapknockoff-master/R/generate_statistics.R")
source("../externalCodes/cheapknockoff-master/R/gz.R")

cheapko_filter <- function(X, Xk.list, y, q){
  
  if("list" %in% class(Xk.list)){
    Xk = do.call("cbind", Xk.list) 
    M.ko = length(Xk.list)
  } 
  else if("matrix" %in% class(Xk.list) | "dataframe" %in% class(Xk.list)){
    Xk = Xk.list 
    M.ko = ncol(Xk)/ncol(X)
  } 
  omega <- rep(M.ko+1, ncol(X))
  
  mko.stat <- stat_glmnet_coef(X, Xk, y, omega = omega)
  mko.result <- filter_gz(mko.stat$kappa, mko.stat$tau, fdr = q, n_knockoff = M.ko)
  mko.result$S <- mko.result$selected
  return(mko.result)
}
