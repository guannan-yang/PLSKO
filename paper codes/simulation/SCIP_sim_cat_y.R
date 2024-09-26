require(glmnet)
require(graphics)
library(MASS)
library(Matrix)
library(knockoff)


library(tidyverse)
library(mixOmics)
library(CVglasso)
library(glmnet)
library(reshape2)
library(irlba)
library(kernlab) # MMD measurement

#library(clustermq)
#library(PerMallows)
library(rARPACK)
library(doParallel)

#  cl <- startMPIcluster(verbose=TRUE) # by default will start one fewer slave
# # # than elements in tasks requested (to account for master)
#  registerDoMPI(cl)

setwd("./plsko_writing/paper codes")

sim_exp_caty <- function(n, X.method = "fac2", p.in.cluster, p.indep = 0, cluster = 1, rho = 0.5, r.in.cluster = 3, fac.load = c(1,2), x.noise = 1,
                         s = 0.1, a = c(3,5), y.dis = "Normal", q = 0.05,
                         plsko.threshold.q = 0.8, full.ncomp = r.in.cluster*cluster-1, plsko.ncomp = r.in.cluster,
                         seed.set = Sys.time(), call.var = "seed.set"){


  source('./RSource/data_sim1.R')
  source('./RSource/DGPX.R')
  source('./RSource/plsko_adv.R')
  source('./RSource/plsko_plsonly.R')
  source('./RSource/KO_perf.R')
  source('./RSource/KOBT/create.pc.knockoff.R')
  source('./RSource/hdko.R')

  set.seed(seed.set)

  GX <- DGPX.gauss(n, X.method = X.method, p.in.cluster = p.in.cluster, p.indep = p.indep, cluster = cluster,
                   rho = rho, r.in.cluster = r.in.cluster, load = fac.load, noise = x.noise, seed = seed.set)
  X <- GX$X
  X <- scale(X, center = F)
  true.cov <- GX$cov
  if(!any(is.na(true.cov))){
    true.cov <- cov2cor(true.cov)
  }

  p <- ncol(X)
  s = ceiling(s*p)
  Gy <- DGPy.AR1(X, s = s, A = a, c = 0, y.dis = "Binary")
  y <- Gy$y

  JS.ko <- knockoffX.sim(X, method ="JS")
  ipad.ko <- knockoffX.sim(X, method = "ipad.pca.center", r = full.ncomp)
  #ginv.ko <- hd.ko(X, fit.method = "ginv", threshold.q = plsko.threshold.q, r.fit = plsko.ncomp)
  pca.ko <- hd.ko(X, fit.method = "pca", threshold.q = plsko.threshold.q, r.fit = plsko.ncomp)
  plsko <- hd.ko(X, fit.method = "pls", threshold.q = plsko.threshold.q, r.fit = plsko.ncomp)
  perm <- knockoffX.sim(X, method = "perm")

  Xk <- list(JS = JS.ko$ko,
             IPAD = ipad.ko$ko,
             #Ginv = ginv.ko,
             PC = pca.ko,
             PLSKO = plsko,
             Permuted = perm$ko
  )

  knockx_set <- names(Xk)

  if(!any(is.na(true.cov))){
    mu <- colMeans(X)
    oracle.ko <- create.gaussian(X, mu, true.cov)
    Xk[[length(Xk)+1]] <- oracle.ko

    names(Xk) <- c(knockx_set, "oracle")
  }

  S <- lapply(Xk, function(ko){
    warning("1")
    selected <- ko.filter(X = X, Xk = ko, y = y, q = q, method = "lasso.logistic")
  })


  fdr <- lapply(S, function(s){fdr = fdr(s$S, Gy$Beta)})
  names(fdr) <- names(Xk)
  fdr.plus <- lapply(S, function(s){fdr.plus = fdr(s$S.plus, Gy$Beta)})
  names(fdr.plus) <- names(Xk)

  power <- lapply(S, function(s){pow = pow(s$S, Gy$Beta)})
  names(power) <- names(Xk)
  power.plus <- lapply(S, function(s){pow = pow(s$S.plus, Gy$Beta)})
  names(power.plus) <- names(Xk)


  # call.var
  call <- mget(call.var)
  call <- paste(unlist(call), collapse = ",")

  if(!any(is.na(true.cov))){
    result <- data.frame(ko.method = c(knockx_set, "oracle"),
                         var = unlist(call),
                         fdp = unlist(fdr),
                         tpp = unlist(power),
                         fdp.plus = unlist(fdr.plus),
                         tpp.plus = unlist(power.plus))
  }
  else{
    result <- data.frame(ko.method = knockx_set,
                         var = unlist(call),
                         fdp = unlist(fdr),
                         tpp = unlist(power),
                         fdp.plus = unlist(fdr.plus),
                         tpp.plus = unlist(power.plus))
  }

  BH.selected <- S_BH(X, y, q = q)

  BH.result <- data.frame(ko.method = "BH",
                          var = unlist(call),
                          fdp = fdr(BH.selected, Gy$Beta),
                          tpp = pow(BH.selected, Gy$Beta),
                          fdp.plus = NA,
                          tpp.plus = NA)
  result <- rbind(result, BH.result)

  return(result)
}



#cl <- makeCluster(detectCores())
cl <- makeCluster(11)
registerDoParallel(cl)
it = 100

result_caty_n_fac2 <- foreach(i = rep(1:it,9), N = rep(seq(200, 1000, 100), each = it),
                              .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                res <- sim_exp_caty(n = N, X.method = "fac2", p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                              }
result_caty_n_fac2_df <- do.call("rbind",result_caty_n_fac2)
saveRDS(result_caty_n_fac2_df, "./linear_y_sim/scip_caty_n_fac2.rds")


result_caty_n_equi5 <- foreach(i = rep(1:it,9), N = rep(seq(200, 1000, 100), each = it),
                               .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                 res <- sim_exp_caty(n = N, X.method = "equi", p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                               }
result_caty_n_equi5_df <- do.call("rbind",result_caty_n_equi5)
saveRDS(result_caty_n_equi5_df, "./linear_y_sim/scip_caty_n_equi5.rds")
