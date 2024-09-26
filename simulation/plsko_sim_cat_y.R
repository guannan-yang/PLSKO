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
library(rARPACK)
library(doParallel)

setwd("./plsko_writing/paper codes")

sim_exp_caty_seqko <- function(n, X.method = "fac2", p.in.cluster, p.indep = 0, cluster = 1, rho = 0.5, r.in.cluster = 3, fac.load = c(1,2), x.noise = 1,
                               s = 0.1, a = c(3,5), y.dis = "Normal", q = 0.05,
                               plsko.threshold.q = 0.8, full.ncomp = r.in.cluster*cluster-1, plsko.ncomp = r.in.cluster,
                               seed.set = Sys.time(), call.var = "seed.set"){

  source('./RSource/data_sim1.R')
  source('./RSource/DGPX.R')
  source('./RSource/plsko_plsonly.R')
  source("./RSource/seqknockoff/internal.R")
  source("./RSource/seqknockoff/knockoff_filters.R")
  source("./RSource/seqknockoff/simknockoffs.R")
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

  seq.ko <- knockoffX.sim(X, method = "seqknockoffs", r = plsko.ncomp)
  pca.ko.orig <- knockoffX.sim(X, method = "PCKO", r = plsko.ncomp, threshold.q = 0)
  pca.ko.orig.full <- knockoffX.sim(X, method = "PCKO", r = full.ncomp, threshold.q = 0)

  # if("fac1" %in% X.method | "fac2" %in% X.method) linear.thre = plsko.ncomp
  # else linear.thre = nrow(X)/2

  plsko.plsonly <- knockoffX.sim(X, method = "PLSKO.only", threshold.q = plsko.threshold.q, r = plsko.ncomp)
  plsko.plsonly.full <- knockoffX.sim(X, method = "PLSKO.only", threshold.q = 0, r = plsko.ncomp)
  plsko.plsonly.full.sparse <- knockoffX.sim(X, method = "PLSKO.only", threshold.abs = 0, r = plsko.ncomp, sparsity = 1/cluster)

  Xk <- list(seqko = seq.ko$ko,
             pcko.orig = pca.ko.orig$ko,
             pcko.orig.full = pca.ko.orig.full$ko,
             plsonly = plsko.plsonly$ko,
             plsonly.full = plsko.plsonly.full$ko,
             plsonly.full.sparse = plsko.plsonly.full.sparse$ko
  )
  knockx_set <- names(Xk)

  run.time <- list(seqko = seq.ko$run.time,
                   pcko.orig = pca.ko.orig$run.time,
                   pcko.orig.full = pca.ko.orig.full$run.time,
                   plsonly = plsko.plsonly$run.time,
                   plsonly.full = plsko.plsonly.full$run.time,
                   plsonly.full.sparse = plsko.plsonly.full.sparse$run.time
  )

  S <- lapply(Xk, function(ko){
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

  result <- data.frame(ko.method = knockx_set,
                       var = unlist(call),
                       fdp = unlist(fdr),
                       tpp = unlist(power),
                       fdp.plus = unlist(fdr.plus),
                       tpp.plus = unlist(power.plus),
                       runtime = unlist(run.time))


  return(result)
}



#cl <- makeCluster(detectCores())
cl <- makeCluster(11)
registerDoParallel(cl)
it = 100

plsko_caty_n_fac2 <- foreach(i = rep(1:it,9), N = rep(seq(200, 1000, 100), each = it),
                             .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                               res <- sim_exp_caty_seqko(n = N, X.method = "fac2", p.in.cluster = 100, cluster = 5, seed.set = i,
                                                         x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                             }
plsko_caty_n_fac2_df <- do.call("rbind",plsko_caty_n_fac2)
saveRDS(plsko_caty_n_fac2_df, "./linear_y_sim/plsko_caty_n_fac2.rds")


plsko_caty_n_equi5 <- foreach(i = rep(1:it,9), N = rep(seq(200, 1000, 100), each = it),
                              .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                res <- sim_exp_caty_seqko(n = N, X.method = "equi", p.in.cluster = 100, cluster = 5, seed.set = i,
                                                          x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                              }
plsko_caty_n_equi5_df <- do.call("rbind",plsko_caty_n_equi5)
saveRDS(plsko_caty_n_equi5_df, "./linear_y_sim/plsko_caty_n_equi5.rds")
