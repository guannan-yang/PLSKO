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

sim_exp_cluster <- function(n, X.method = "fac2", p.in.cluster, p.indep = 0, cluster = 1, rho = 0.5, r.in.cluster = 3,  x.noise = 1,
                            s = 0.1, a = c(3,5), y.dis = "Normal", q = 0.05, y.noise = 0,
                            plsko.threshold.q = 0.8, full.ncomp = r.in.cluster*cluster-1, ko.diagnosis = F, plsko.ncomp = r.in.cluster,
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
                   rho = rho, r.in.cluster = r.in.cluster,  noise = x.noise, seed = seed.set)
  X <- GX$X
  X <- scale(X, center = F)
  true.cov <- GX$cov
  true.cov <- cov2cor(true.cov)

  p <- ncol(X)
  s = ceiling(s*p)
  Gy <- DGPy.AR1(X, s = s, A = a, c = y.noise, y.dis = "Normal")
  y <- Gy$y


  seq.ko <- knockoffX.sim(X, method = "seqknockoffs", ko.diagnosis = ko.diagnosis, r = plsko.ncomp)
  pca.ko.orig <- knockoffX.sim(X, method = "PCKO", ko.diagnosis = ko.diagnosis, r = plsko.ncomp, threshold.q = 0)
  pca.ko.orig.full <- knockoffX.sim(X, method = "PCKO", ko.diagnosis = ko.diagnosis, r = plsko.ncomp, threshold.q = 0)

  # if("fac1" %in% X.method | "fac2" %in% X.method) linear.thre = plsko.ncomp
  # else linear.thre = nrow(X)/2

  plsko.plsonly <- knockoffX.sim(X, method = "PLSKO.only", threshold.q = plsko.threshold.q, ko.diagnosis = ko.diagnosis, r = plsko.ncomp)
  plsko.plsonly.full <- knockoffX.sim(X, method = "PLSKO.only", threshold.q = 0, r = plsko.ncomp, ko.diagnosis = ko.diagnosis)
  plsko.plsonly.full.sparse <- knockoffX.sim(X, method = "PLSKO.only", threshold.abs = 0, r = plsko.ncomp, sparsity = 1/cluster, ko.diagnosis = ko.diagnosis)

  Xk <- list(seqko = seq.ko$ko,
             pcko.orig = pca.ko.orig$ko,
             pcko.orig.full = pca.ko.orig.full$ko,
             plsonly = plsko.plsonly$ko,
             plsonly.full = plsko.plsonly.full$ko,
             plsonly.full.sparse = plsko.plsonly.full.sparse$ko
  )
  knockx_set <- names(Xk)

  #diagnosis
  if(ko.diagnosis){
    diagnosis <- list(seqko = seq.ko$diagnosis,
                      pcko.orig = pca.ko.orig$diagnosis,
                      pcko.orig.full = pca.ko.orig.full$diagnosis,
                      plsonly = plsko.plsonly$diagnosis,
                      plsonly.full = plsko.plsonly.full$diagnosis,
                      plsonly.full.sparse = plsko.plsonly.full.sparse$ko
    )
  }

  run.time <- list(seqko = seq.ko$run.time,
                   pcko.orig = pca.ko.orig$run.time,
                   pcko.orig.full = pca.ko.orig.full$run.time,
                   plsonly = plsko.plsonly$run.time,
                   plsonly.full = plsko.plsonly.full$run.time,
                   plsonly.full.sparse = plsko.plsonly.full.sparse$run.time
  )


  if(ko.diagnosis){
    diagnosis.df <- do.call("rbind", diagnosis)
  }


  S <- lapply(Xk, function(ko){
    selected <- ko.filter(X = X, Xk = ko, y = y, q = q, method = "lasso.lcd")
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

  result <- data.frame(ko.method = knockx_set,
                       X.method = X.method,
                       var = unlist(call),
                       fdp = unlist(fdr),
                       tpp = unlist(power),
                       fdp.plus = unlist(fdr.plus),
                       tpp.plus = unlist(power.plus),
                       run.time = unlist(run.time))
  #}

  if(ko.diagnosis){
    result <- cbind(result, diagnosis.df)
  }

  return(result)
}

#cl <- makeCluster(detectCores())
cl <- makeCluster(11)
registerDoParallel(cl)
it = 100

## series
#N
result_plsko_n <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                          .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                            res <- sim_exp_cluster(n = N, X.method = "fac2", p.in.cluster = 100, cluster = 5, s = 0.1, seed.set = i, call.var = "n")
                          }
result_plsko_n_df <- do.call("rbind", result_plsko_n)
saveRDS(result_plsko_n_df, "./linear_y_sim/result_plsko_n.rds")

#ynoise
result_plsko_noisey_fac2<- foreach(i = rep(1:it,5), ynoise = rep(1:5, each = it),
                                   .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                     res <- sim_exp_cluster(n = 500, X.method = "fac2", y.noise = ynoise, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                            plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "y.noise")
                                   }
result_plsko_noisey_fac2_df <- do.call("rbind",result_plsko_noisey_fac2)
saveRDS(result_plsko_noisey_fac2_df, "./linear_y_sim/plsko_noisey_fac2.rds")


# sparsity
result_plsko_s_fac2 <- foreach(i = rep(1:it,4), sparsity = rep(seq(0.1,0.4,0.1), each = it),
                               .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                 res <- sim_exp_cluster(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                        s = sparsity, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "s")
                               }
result_plsko_s_fac2_df <- do.call("rbind",result_plsko_s_fac2)
saveRDS(result_plsko_s_fac2_df, "./linear_y_sim/plsko_s_fac2.rds")


result_plsko_n_equi5 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_cluster(n = N, X.method = "equi", rho = 0.5, p.in.cluster = 100, cluster = 5,
                                                         x.noise = 0,s = 0.1, seed.set = i, call.var = "n")
                                }
result_plsko_n_equi5_df <- do.call("rbind", result_plsko_n_equi5)
saveRDS(result_plsko_n_equi5_df, "./linear_y_sim/result_plsko_n_equi5.rds")

result_plsko_n_ar1 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                              .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                res <- sim_exp_cluster(n = N, X.method = "AR1", rho = 0.5, p.in.cluster = 100, cluster = 5,
                                                       x.noise = 0,s = 0.1, seed.set = i, call.var = "n")
                              }
result_plsko_n_ar1_df <- do.call("rbind", result_plsko_n_ar1)
saveRDS(result_plsko_n_ar1_df, "./linear_y_sim/result_plsko_n_ar1.rds")


result_plsko_n_fac2quad <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                   .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                     res <- sim_exp_cluster(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                            x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 6, call.var = c("n", "seed.set"))
                                   }
result_plsko_n_fac2quad_df <- do.call("rbind",result_plsko_n_fac2quad)
saveRDS(result_plsko_n_fac2quad_df, "./linear_y_sim/result_plsko_n_fac2quad.rds")

result_plsko_n_fac2quad2 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                    .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                      res <- sim_exp_cluster(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                             x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 9, call.var = c("n", "seed.set"))
                                    }
result_plsko_n_fac2quad2_df <- do.call("rbind",result_plsko_n_fac2quad2)
saveRDS(result_plsko_n_fac2quad2_df, "./linear_y_sim/result_plsko_n_fac2quad2.rds")

result_plsko_n_fac2quad3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                    .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                      res <- sim_exp_cluster(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                             x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n", "seed.set"))
                                    }
result_plsko_n_fac2quad3_df <- do.call("rbind",result_plsko_n_fac2quad3)
saveRDS(result_plsko_n_fac2quad3_df, "./linear_y_sim/result_plsko_n_fac2quad3.rds")
