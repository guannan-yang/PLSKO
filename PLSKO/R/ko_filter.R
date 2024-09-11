
ko.filter <- function(X, Xk, y, q, method = "lasso.lcd", offset = 0, ...){
  n = nrow(X)
  p = ncol(X)
  X.names = names(X)

  if(method == "lasso.lcd"){
    W <- stat.lasso_coefdiff(X, Xk, y, ...)
  }
  else if(method == "lasso.max.lambda"){
    W <- stat.lasso_lambdadiff(X, Xk, y, ...)
  }
  # else if(method == "pls.regression.lcd"){
  #   W <- W_spls_lcd(X, Xk, y, lcd_comp = 1)
  # }
  #
  # else if(method == "plsda.lcd"){ # exclude when continuous y design
  #   W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1)
  # }
  # else if(method == "spls.regression.lcd"){
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_spls_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }
  # else if(method == "splsda.lcd"){ # exclude when continuous y design
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_splsda_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }

  else if(method == "lasso.logistic"){ # exclude when continuous y design
    W <- stat.lasso_coefdiff_bin(X, Xk, y, ...)
  }

  else if(method == "RF"){
    W <- stat.random_forest(X, Xk, y)
  }

  # find the knockoff threshold T
  t = sort(c(0, abs(W)))

  ratio = c(rep(0, p))
  ratio.plus = c(rep(0, p))

  for (j in 1:p) {
    if(offset==0){
    ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }
    else if(offset == 1){
    ratio[j] = (1+sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }
  }

  id = which(ratio <= q)[1]

  if(length(id) == 0){
    T = Inf
  } else {
    T = t[id]
  }

  # set of discovered variables
  S = integer0_test(which(W >= T))

  return(list(S = S, W = W))
}

ko.withW <- function()


  
integer0_test <- function(which){
  if(identical(which, integer(0))){
    return(NULL)
  }
  else {return(which)}
}
