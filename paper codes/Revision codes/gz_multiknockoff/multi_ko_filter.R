# gogohugo cheapknockoff filtering ======
library(CVXR)
# 
# source("../externalCodes/cheapknockoff-master/R/generate_statistics.R")
# source("../externalCodes/cheapknockoff-master/R/gz.R")

cheapko_filter <- function(X, Xk.list, y, q, family = "gaussian", offset = 1){
  
  if("list" %in% class(Xk.list)){
    Xk = do.call("cbind", Xk.list) 
    M.ko = length(Xk.list)
  } 
  else if("matrix" %in% class(Xk.list) | "dataframe" %in% class(Xk.list)){
    Xk = Xk.list 
    M.ko = ncol(Xk)/ncol(X)
  } 
  omega <- rep(M.ko+1, ncol(X))
  
  mko.stat <- stat_glmnet_coef(X, Xk, y, omega = omega, family = family)
  mko.result <- filter_gz(mko.stat$kappa, mko.stat$tau, fdr = q, n_knockoff = M.ko, offset = offset)
  #mko.result$S <- mko.result$selected
  return(mko.result)
}
