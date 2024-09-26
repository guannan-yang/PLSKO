library(JuliaCall)
#julia <- julia_setup()
#julia_install_package_if_needed("Knockoffs")
julia_library("Knockoffs")

generate_ko_with_julia <- function(X, Sigma = NULL, mu, method, m = 1, windowsize = 100) {
  start.time <- Sys.time()
  # pass variables from R to Julia
  julia_assign("X", X)
  julia_assign("Sigma", Sigma)
  julia_assign("mu", mu)
  julia_assign("method", method)
  julia_assign("m", m)
  # with Julia, solve group knockoffs problem
  #julia_command("result = modelX_gaussian_knockoffs(X, Symbol(method), mu, Sigma, m=Int(m), verbose=true)")
  
  if(is.null(Sigma)|any(is.na(Sigma))){ # using shrinked covariance by ledoit-wolf
    if(method == "mvr"){
      julia_command("result = modelX_gaussian_knockoffs(X, Symbol(:mvr))")
    }
    else if(method == "maxent"){
      julia_command("result = modelX_gaussian_knockoffs(X, Symbol(:maxent))")
    }
    else if(method == "approx"){
      julia_assign("windowsize", windowsize)
      julia_command("result = approx_modelX_gaussian_knockoffs(X, Symbol(:maxent), windowsize=Int(windowsize))")
    }
  }
  else{ # using oracle true covariance 
     if(method == "mvr"){
    julia_command("result = modelX_gaussian_knockoffs(X, Symbol(:mvr), mu, Sigma, m=Int(m), verbose=true)")
  }
  else if(method == "maxent"){
    julia_command("result = modelX_gaussian_knockoffs(X, Symbol(:maxent), mu, Sigma, m=Int(m), verbose=true)")
  } 
  }


  # pull variables from Julia back to R
  #Xko <- julia_eval("result.X̃ ") # the knockoff
  Xko <- julia_eval("result.Xko") # the knockoff
  end.time <- Sys.time()
  run.time <- difftime(end.time, start.time, units = "secs")
  
  result <- list("Xko"=Xko, 
                 "run.time"=run.time)
  return (result)
}

# setwd("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2")
# 
# source('./RSource/data_sim1.R')
# source('./RSource/DGPX.R')
# 
# test.x <- DGPX.gauss(500,100,X.method = "equi", cluster = 1, noise = 0.5, seed = 1)
# test.x.mat <- scale(test.x$X, center = F)
# test.y <- DGPy.AR1(test.x.mat, s = ceiling(0.1*ncol(test.x.mat)), c = 0, A = c(3,5))
# 
# test <- generate_ko_with_julia(X = test.x.mat, Sigma = cov2cor(test.x$cov), mu = rep(0, 100), method ="mvr",m = 1)
