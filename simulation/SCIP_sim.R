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

#  cl <- startMPIcluster(verbose=TRUE) # by default will start one fewer slave
# # # than elements in tasks requested (to account for master)
#  registerDoMPI(cl)

setwd("./plsko_writing/paper codes")

sim_exp_hd <- function(n, X.method = "fac2", p.in.cluster, p.indep = 0, cluster = 1, rho = 0.5, r.in.cluster = 3, fac.load = c(1,2), x.noise = 1,
                       s = 0.1, a = c(3,5), y.dis = "Normal", q = 0.05, y.noise = 0,
                       plsko.threshold.q = 0.8, full.ncomp = r.in.cluster*cluster-1, plsko.ncomp = r.in.cluster,
                       seed.set = Sys.time(), call.var = "seed.set"){


  source('./RSource/data_sim1.R')
  source('./RSource/DGPX.R')
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
  Gy <- DGPy.AR1(X, s = s, A = a, c = y.noise, y.dis = "Normal")
  y <- Gy$y

  JS.ko <- knockoffX.sim(X, method ="JS")
  ipad.ko <- knockoffX.sim(X, method = "ipad.pca.center", r = full.ncomp)
  ipad.ko.full <- knockoffX.sim(X, method = "ipad.pca.center", r = plsko.ncomp * cluster)

  #ginv.ko <- hd.ko(X, fit.method = "ginv", threshold.q = plsko.threshold.q, r.fit = plsko.ncomp)
  pca.ko <- hd.ko(X, fit.method = "pca", threshold.q = plsko.threshold.q, r = plsko.ncomp)
  plsko <- hd.ko(X, fit.method = "pls", threshold.q = plsko.threshold.q, r = plsko.ncomp)
  #perm <- knockoffX.sim(X, method = "perm")

  Xk <- list(JS = JS.ko$ko,
             IPAD = ipad.ko$ko,
             #Ginv = ginv.ko,
             PC = pca.ko,
             PLSKO = plsko,
             IPAD.full = ipad.ko.full$ko
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
    selected <- ko.filter(X = X, Xk = ko, y = y, q = q, method = "lasso.lcd")
  })
  #names(S) <- c(knockx_set, "oracle")


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
                         X.method = X.method,
                         var = unlist(call),
                         fdp = unlist(fdr),
                         tpp = unlist(power),
                         fdp.plus = unlist(fdr.plus),
                         tpp.plus = unlist(power.plus))
  }
  else{
    result <- data.frame(ko.method = knockx_set,
                         X.method = X.method,
                         var = unlist(call),
                         fdp = unlist(fdr),
                         tpp = unlist(power),
                         fdp.plus = unlist(fdr.plus),
                         tpp.plus = unlist(power.plus))
  }

  return(result)
}

#cl <- makeCluster(detectCores())
cl <- makeCluster(11)
registerDoParallel(cl)
it = 100




### series
result_hd_n_fac2 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                              res <- sim_exp_hd(n = N, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                            }
result_hd_n_fac2_df <- do.call("rbind",result_hd_n_fac2)
saveRDS(result_hd_n_fac2_df, "./linear_y_sim/oracle/scip_n_fac2.rds")


result_hd_noisey_fac2<- foreach(i = rep(1:it,5), ynoise = rep(1:5, each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = ynoise, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "y.noise")
                                }
result_hd_noisey_fac2_df <- do.call("rbind",result_hd_noisey_fac2)
saveRDS(result_hd_noisey_fac2_df, "./linear_y_sim/oracle/scip_noisey_fac2.rds")


result_hd_ncluster_lown250<- foreach(i = rep(1:it,5), ncluster = rep(c(1,2,4,5,10), each = it), pincluster = rep(c(500,250,125,100,50), each = it),
                                     thres.q = rep(c(0,0.5,0.75,0.8,0.9), each = it),
                                     .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                       res <- sim_exp_hd(n = 250, X.method = "fac2",rho = 0.7, y.noise = 0, p.in.cluster = pincluster, cluster = ncluster, seed.set = i,
                                                         plsko.threshold.q = thres.q, plsko.ncomp = 3, call.var = "cluster")
                                     }
result_hd_ncluster_lown250_df <- do.call("rbind",result_hd_ncluster_lown250)
saveRDS(result_hd_ncluster_lown250_df, "./linear_y_sim/oracle/scip_ncluster_lown250.rds")

### sparsity when N = 500
result_hd_s_fac2<- foreach(i = rep(1:it,4), sparsity = rep(seq(0.1,0.4,0.1), each = it),
                           .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                             res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                               s = sparsity, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "s")
                           }
result_hd_s_fac2_df <- do.call("rbind",result_hd_s_fac2)
saveRDS(result_hd_s_fac2_df, "./linear_y_sim/oracle/scip_s.rds")

result_hd_s_fac2_lown350<- foreach(i = rep(1:it,4), sparsity = rep(seq(0.1,0.4,0.1), each = it),
                                   .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                     res <- sim_exp_hd(n = 350, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                       s = sparsity, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "s")
                                   }
result_hd_s_fac2_lown350_df <- do.call("rbind",result_hd_s_fac2_lown350)
saveRDS(result_hd_s_fac2_lown350_df, "./linear_y_sim/scip_s_lown350.rds")




#### fac1
result_hd_n_fac1 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                              res <- sim_exp_hd(n = N, X.method = "fac1", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                            }
result_hd_n_fac1_df <- do.call("rbind",result_hd_n_fac1)
saveRDS(result_hd_n_fac1_df, "./linear_y_sim/oracle/scip_n_fac1.rds")


result_hd_noisey_fac1<- foreach(i = rep(1:it,5), ynoise = rep(1:5, each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = 500, X.method = "fac1", y.noise = ynoise, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "y.noise")
                                }
result_hd_noisey_fac1_df <- do.call("rbind",result_hd_noisey_fac1)
saveRDS(result_hd_noisey_fac1_df, "./linear_y_sim/oracle/scip_noisey_fac1.rds")

result_hd_s_fac1<- foreach(i = rep(1:it,4), sparsity = rep(seq(0.1,0.4,0.1), each = it),
                           .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                             res <- sim_exp_hd(n = 500, X.method = "fac1", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                               s = sparsity, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "s")
                           }
result_hd_s_fac1_df <- do.call("rbind",result_hd_s_fac1)
saveRDS(result_hd_s_fac1_df, "./linear_y_sim/oracle/scip_s_fac1.rds")


## fac2 with different 'sharper' load with b:a = 1.2
result_hd_n_fac2_alterload <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                                      .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                        res <- sim_exp_hd(n = N, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                          x.noise = 1, fac.load = c(1, 1.2), plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                                      }
result_hd_n_fac2_alterload_df <- do.call("rbind",result_hd_n_fac2_alterload)
saveRDS(result_hd_n_fac2_alterload_df, "./linear_y_sim/oracle/scip_n_fac2_alterload.rds")

## robustness test
# wrong number of component
result_hd_comp_fac2<- foreach(i = rep(1:it,4), ko.ncomp = rep(2:5, each = it),
                              .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                  plsko.threshold.q = 0.8, plsko.ncomp = ko.ncomp, call.var = "plsko.ncomp")
                              }
result_hd_comp_fac2_df <- do.call("rbind",result_hd_comp_fac2)
saveRDS(result_hd_comp_fac2_df, "./linear_y_sim/scip_comp.rds")

# wrong threshold
result_hd_thres_fac2<- foreach(i = rep(1:it,5), ko.thres = rep(seq(0.5,0.9,0.1), each = it),
                               .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                 res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                   plsko.threshold.q = ko.thres, plsko.ncomp = 3, call.var = "plsko.ncomp") # call.var set wrong..
                               }
result_hd_thres_fac2_df <- do.call("rbind",result_hd_thres_fac2)
#result_hd_noisey_fac2_df$var <- rep(seq(0.5,0.9,0.1), each = it)
saveRDS(result_hd_thres_fac2_df, "./linear_y_sim/scip_thres.rds")

# wrong threshold and component
result_hd_thres_ncomp_fac2<- foreach(i = rep(1:it,4), ko.thres = rep(seq(0.6,0.9,0.1), it), ko.ncomp = rep(2:5, each = it),
                                     .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                       res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                         plsko.threshold.q = ko.thres, plsko.ncomp = ko.ncomp, call.var = c("plsko.ncomp", "plsko.threshold.q"))
                                     }
result_hd_thres_ncomp_fac2_df <- do.call("rbind",result_hd_thres_ncomp_fac2)
saveRDS(result_hd_thres_ncomp_fac2_df, "./linear_y_sim/scip_thres_ncomp.rds")


#### robust more extreme
# wrong number of component
result_hd_comp_fac2_2<- foreach(i = rep(1:it,4), ko.ncomp = rep(c(6:9), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    plsko.threshold.q = 0.8, plsko.ncomp = ko.ncomp, call.var = "plsko.ncomp")
                                }
result_hd_comp_fac2_2_df <- do.call("rbind",result_hd_comp_fac2_2)
saveRDS(result_hd_comp_fac2_2_df, "./linear_y_sim/scip_comp_2.rds")

# wrong threshold
result_hd_thres_fac2_2<- foreach(i = rep(1:it,5), ko.thres = rep(seq(0,0.4,0.1), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                   res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     plsko.threshold.q = ko.thres, plsko.ncomp = 3, call.var = "plsko.threshold.q") # call.var set wrong..
                                 }
result_hd_thres_fac2_2_df <- do.call("rbind",result_hd_thres_fac2_2)
saveRDS(result_hd_thres_fac2_2_df, "./linear_y_sim/scip_thres_2.rds")

### Xnoise lower
result_hd_xnoise_fac2<- foreach(i = rep(1:it,4), xnoise = rep(c(0.1,0.5,1,2), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    x.noise = xnoise, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "x.noise")
                                }
result_hd_xnoise_fac2_df <- do.call("rbind",result_hd_xnoise_fac2)
saveRDS(result_hd_xnoise_fac2_df, "./linear_y_sim/scip_xnoise.rds")

### rank
result_hd_r_fac2<- foreach(i = rep(1:it,5), X.r = rep(1:5, each = it),
                           .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                             res <- sim_exp_hd(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                               r.in.cluster = X.r, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "r.in.cluster")
                           }
result_hd_r_fac2_df <- do.call("rbind",result_hd_r_fac2)
saveRDS(result_hd_r_fac2_df, "./linear_y_sim/scip_r.rds")


result_hd_rho_lown250<- foreach(i = rep(1:it,5), rho = rep(seq(0.1,0.9,0.2), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = 250, X.method = "equi", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    x.noise = 0, rho = rho, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "rho")
                                }
result_hd_rho_lown250_df <- do.call("rbind",result_hd_rho_lown250)
saveRDS(result_hd_rho_lown250_df, "./linear_y_sim/scip_rho_lown250.rds")

### AR1 rho
result_hd_rho_lown250_ar1 <- foreach(i = rep(1:it,5), rho = rep(seq(0.1,0.9,0.2), each = it),
                                     .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                       res <- sim_exp_hd(n = 250, X.method = "AR1", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                         x.noise = 0, rho = rho, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "rho")
                                     }
result_hd_rho_lown250_ar1_df <- do.call("rbind",result_hd_rho_lown250_ar1)
saveRDS(result_hd_rho_lown250_ar1_df, "./linear_y_sim/scip_rho_lown250_ar1.rds")

### AR1 N
result_hd_n_ar1 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                           .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                             res <- sim_exp_hd(n = N, X.method = "AR1", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                               x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                           }
result_hd_n_ar1_df <- do.call("rbind",result_hd_n_ar1)
saveRDS(result_hd_n_ar1_df, "./linear_y_sim/scip_n_ar1.rds")

### equi N
result_hd_n_equi <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                              res <- sim_exp_hd(n = N, X.method = "equi", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                            }
result_hd_n_equi_df <- do.call("rbind",result_hd_n_equi)
saveRDS(result_hd_n_equi_df, "./linear_y_sim/oracle/scip_n_equi.rds")


# Quadratic X, N series
result_hd_n_quad <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                              tryCatch({res <- sim_exp_hd(n = N, X.method = "quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                          x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 6, call.var = c("n", "seed.set"))},
                                       error = function(e){cat("ERROR :",conditionMessage(e)," ",i," ",N, "\n")})

                            }
result_hd_n_quad_df <- do.call("rbind",result_hd_n_quad)
saveRDS(result_hd_n_quad_df, "./linear_y_sim/scip_n_quad.rds")

result_hd_n_fac2quad <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 6, call.var = c("n", "seed.set"))
                                }
result_hd_n_fac2quad_df <- do.call("rbind",result_hd_n_fac2quad)
saveRDS(result_hd_n_fac2quad_df, "./linear_y_sim/scip_n_fac2quad.rds")

result_hd_n_fac2quad2 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                   res <- sim_exp_hd(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 9, call.var = c("n", "seed.set"))
                                 }
result_hd_n_fac2quad2_df <- do.call("rbind",result_hd_n_fac2quad2)
saveRDS(result_hd_n_fac2quad2_df, "./linear_y_sim/scip_n_fac2quad2.rds")

result_hd_n_fac2quad3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                   res <- sim_exp_hd(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n", "seed.set"))
                                 }
result_hd_n_fac2quad3_df <- do.call("rbind",result_hd_n_fac2quad3)
saveRDS(result_hd_n_fac2quad3_df, "./linear_y_sim/scip_n_fac2quad3.rds")

# fac2-sin ======
result_hd_n_fac2sin3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                  res <- sim_exp_hd(n = N, X.method = "fac2-sin", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                    x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n", "seed.set"))
                                }
result_hd_n_fac2sin3_df <- do.call("rbind",result_hd_n_fac2sin3)
saveRDS(result_hd_n_fac2sin3_df, "./linear_y_sim/scip_n_fac2sin3.rds")

#mixed X
result_hd_n_mix <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                           .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                             res <- sim_exp_hd(n = N, X.method = c("fac1","fac2","equi","ar1","quad"), y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                               x.noise = 0.5, rho = 0.8, plsko.threshold.q = 0.8, plsko.ncomp = 3, call.var = "n")
                           }
result_hd_n_mix_df <- do.call("rbind",result_hd_n_mix)
saveRDS(result_hd_n_mix_df, "./linear_y_sim/scip_n_mix.rds")

# fac2 interaction ====
result_hd_n_fac2inter2 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                    method <- "fac2-inter"
                                    res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 9, call.var = c("n", "seed.set"))
                                  }
result_hd_n_fac2inter2_df <- do.call("rbind",result_hd_n_fac2inter2)
saveRDS(result_hd_n_fac2inter2_df, "./linear_y_sim/scip_n_fac2inter2.rds")

result_hd_n_fac2inter3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                    method <- "fac2-inter"
                                    res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 1, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n", "seed.set"))
                                  }
result_hd_n_fac2inter3_df <- do.call("rbind",result_hd_n_fac2inter3)
saveRDS(result_hd_n_fac2inter3_df, "./linear_y_sim/scip_n_fac2inter3.rds")

# equi-quad =====
result_hd_n_equiquad2 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                   method <- "equi-quad"
                                   res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 9, call.var = c("n", "seed.set"))
                                 }
result_hd_n_equiquad2_df <- do.call("rbind",result_hd_n_equiquad2)
saveRDS(result_hd_n_equiquad2_df, "./linear_y_sim/scip_n_equiquad2.rds")

result_hd_n_equiquad3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                   method <- "equi-quad"
                                   res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                     x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n", "seed.set"))
                                 }
result_hd_n_equiquad3_df <- do.call("rbind",result_hd_n_equiquad3)
saveRDS(result_hd_n_equiquad3_df, "./linear_y_sim/scip_n_equiquad3.rds")

# equi-interact =====
result_hd_n_equiinter2 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                    method <- "equi-inter"
                                    res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 9, call.var = c("n"))
                                  }
result_hd_n_equiinter2_df <- do.call("rbind",result_hd_n_equiinter2)
saveRDS(result_hd_n_equiinter2_df, "./linear_y_sim/scip_n_equiinter2.rds")

result_hd_n_equiinter3 <- foreach(i = rep(1:it,8), N = rep(seq(150,500,50), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                    method <- "equi-inter"
                                    res <- sim_exp_hd(n = N, X.method = method, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 0, rho = 0.5, plsko.threshold.q = 0.8, plsko.ncomp = 12, call.var = c("n"))
                                  }
result_hd_n_equiinter3_df <- do.call("rbind",result_hd_n_equiinter3)
saveRDS(result_hd_n_equiinter3_df, "./linear_y_sim/scip_n_equiinter3.rds")
