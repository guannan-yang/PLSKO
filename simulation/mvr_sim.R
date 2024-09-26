setwd("./plsko_writing/paper codes")

source('./RSource/mvr_julia.R')
julia <- julia_setup()

require(glmnet)
require(graphics)
library(MASS)
library(Matrix)
library(knockoff)
library(tidyverse)
library(mixOmics)
library(reshape2)
library(irlba)
library(kernlab) # MMD measurement

library(rARPACK)
library(doParallel)

cl <- makeCluster(11)
registerDoParallel(cl)

source('./RSource/data_sim1.R')
source('./RSource/DGPX.R')

sim_exp_mvr <- function(n, X.method = "fac2", p.in.cluster, p.indep = 0, cluster = 1, rho = 0.5, r.in.cluster = 3, fac.load = c(1,2), x.noise = 1,
                        s = 0.1, a = c(3,5), y.dis = "Normal", q = 0.05, y.noise = 0,
                        w.method = "lasso.lcd",
                        seed.set = Sys.time(), call.var = "seed.set"){

  source('./RSource/data_sim1.R')
  source('./RSource/DGPX.R')
  set.seed(seed.set)

  GX <- DGPX.gauss(n, X.method = X.method, p.in.cluster = p.in.cluster, p.indep = p.indep, cluster = cluster,
                   rho = rho, r.in.cluster = r.in.cluster, load = fac.load, noise = x.noise, seed = seed.set)
  X <- GX$X
  X <- scale(X, center = F)
  true.cov <- GX$cov
  true.cov <- cov2cor(true.cov)
  mu <- colMeans(X)

  p <- ncol(X)
  s = ceiling(s*p)
  Gy <- DGPy.AR1(X, s = s, A = a, c = y.noise, y.dis = y.dis)
  y <- Gy$y

  #knockoff construct
  mvr_oracle <- generate_ko_with_julia(X, true.cov, mu, method = "mvr", m = 1)
  maxent_oracle <- generate_ko_with_julia(X, true.cov, mu, method = "maxent", m = 1)
  mvr_JS <- generate_ko_with_julia(X, Sigma = NULL, mu, method = "mvr", m = 1)
  maxent_JS <- generate_ko_with_julia(X, Sigma = NULL, mu, method = "maxent", m = 1)
  maxent_approx <- generate_ko_with_julia(X, Sigma = NULL, mu, method ="approx", m = 1, windowsize = p.in.cluster)

  Xk <- list(mvr_oracle = mvr_oracle$Xko,
             maxent_oracle = maxent_oracle$Xko,
             mvr_JS = mvr_JS$Xko,
             maxent_JS = maxent_JS$Xko,
             maxent_approx = maxent_approx$Xko)

  knockx.set <- names(Xk)

  # var selection
  S <- lapply(Xk, function(ko){
    selected <- ko.filter(X = X, Xk = ko, y = y, q = q, method = w.method)
  })

  # measurement
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

  result <- data.frame(ko.method = knockx.set,
                       var = unlist(call),
                       fdp = unlist(fdr),
                       tpp = unlist(power),
                       fdp.plus = unlist(fdr.plus),
                       tpp.plus = unlist(power.plus))

  return(result)
}

##### test 1. N series########
it = 100

result_mvr_n_fac2 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                             .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                               source('./RSource/mvr_julia.R')
                               res <- sim_exp_mvr(n = N, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                  x.noise = 1, rho = 0.5, call.var = "n")
                             }
result_mvr_n_fac2_df <- do.call("rbind",result_mvr_n_fac2)
saveRDS(result_mvr_n_fac2_df, "./linear_y_sim/mvr_n_fac2.rds")
result_mvr_n_fac2_df <- readRDS("./linear_y_sim/mvr_n_fac2.rds")

###### test 2. sparsitiy#######

result_mvr_s_fac2<- foreach(i = rep(1:it,4), sparsity = rep(seq(0.1,0.4,0.1), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                              source('./RSource/mvr_julia.R')

                              tryCatch({
                                res <- sim_exp_mvr(n = 500, X.method = "fac2", x.noise = 1, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                   s = sparsity, call.var = "s")

                              }, error = function(e){cat("ERROR :",conditionMessage(e),"\n","i=",i)})
                            }
result_mvr_s_fac2_df <- do.call("rbind",result_mvr_s_fac2)
saveRDS(result_mvr_s_fac2_df, "./linear_y_sim/mvr_s.rds")


###### test 3. noise y #######
result_mvr_noisey_fac2<- foreach(i = rep(1:it,5), ynoise = rep(1:5, each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                   source('./RSource/mvr_julia.R')

                                   tryCatch({
                                     res <- sim_exp_mvr(n = 500, X.method = "fac2", x.noise = 1, y.noise = ynoise, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                        call.var = "y.noise")

                                   }, error = function(e){cat("ERROR :",conditionMessage(e),"\n","i=",i)})
                                 }
result_mvr_noisey_fac2_df <- do.call("rbind",result_mvr_noisey_fac2)
saveRDS(result_mvr_noisey_fac2_df, "./linear_y_sim/mvr_noisey_fac2.rds")


########## test 4: multi gaussian #########

result_mvr_n_equi5 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                              .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                source('./RSource/mvr_julia.R')
                                res <- sim_exp_mvr(n = N, X.method = "equi", rho = 0.5, p.in.cluster = 100, cluster = 5,
                                                   x.noise = 0,s = 0.1, seed.set = i, call.var = "n")
                              }
result_mvr_n_equi5_df <- do.call("rbind", result_mvr_n_equi5)
saveRDS(result_mvr_n_equi5_df, "./linear_y_sim/result_mvr_n_equi5.rds")

result_mvr_n_ar1 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                            .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                              source('./RSource/mvr_julia.R')
                              res <- sim_exp_mvr(n = N, X.method = "AR1", rho = 0.5, p.in.cluster = 100, cluster = 5,
                                                 x.noise = 0,s = 0.1, seed.set = i, call.var = "n")
                            }
result_mvr_n_ar1_df <- do.call("rbind", result_mvr_n_ar1)
saveRDS(result_mvr_n_ar1_df, "./linear_y_sim/result_mvr_n_ar1.rds")

#stopCluster(cl)

########### test 5: categorical y ##########

result_mvr_n_fac2_caty <- foreach(i = rep(1:it,9), N = rep(seq(200, 1000, 100), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                    source('./RSource/mvr_julia.R')
                                    res <- sim_exp_mvr(n = N, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                       y.dis = "Binary", w.method = "lasso.logistic",
                                                       x.noise = 1, rho = 0.5, call.var = "n")
                                  }
result_mvr_n_fac2_caty_df <- do.call("rbind",result_mvr_n_fac2_caty)
saveRDS(result_mvr_n_fac2_caty_df, "./linear_y_sim/mvr_n_fac2_caty.rds")

######### fac2-quad ########

result_mvr_n_fac2_quad <- foreach(i = rep(1:it,11), N = rep(seq(150,500,50), each = it),
                                  .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                    source('./RSource/mvr_julia.R')
                                    res <- sim_exp_mvr(n = N, X.method = "fac2-quad", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                       x.noise = 1, rho = 0.5, call.var = "n")
                                  }
result_mvr_n_fac2_quad_df <- do.call("rbind",result_mvr_n_fac2_quad)
saveRDS(result_mvr_n_fac2_quad_df, "./linear_y_sim/mvr_n_fac2_quad.rds")

## fac2-sin ===
result_mvr_n_fac2_sin <- foreach(i = rep(1:it,11), N = rep(seq(150,500,50), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                   source('./RSource/mvr_julia.R')
                                   res <- sim_exp_mvr(n = N, X.method = "fac2-sin", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 1, rho = 0.5, call.var = "n")
                                 }
result_mvr_n_fac2_sin_df <- do.call("rbind",result_mvr_n_fac2_sin)

result_mvr_rho_lown250<- foreach(i = rep(1:it,5), rho = rep(seq(0.1,0.9,0.2), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                   source('./RSource/mvr_julia.R')
                                   res <- sim_exp_mvr(n = 250, X.method = "equi", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = 0, rho = rho, call.var = "rho")
                                 }
result_mvr_rho_lown250_df <- do.call("rbind",result_mvr_rho_lown250)
saveRDS(result_mvr_rho_lown250_df, "./linear_y_sim/mvr_rho_lown250.rds")

######### fac2-inter ########

result_mvr_n_fac2_inter <- foreach(i = rep(1:it,11), N = rep(seq(150,500,50), each = it),
                                   .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                     source('./RSource/mvr_julia.R')
                                     X.dis <- "fac2-inter"
                                     res <- sim_exp_mvr(n = N, X.method = X.dis, y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                        x.noise = 1, rho = 0.5, call.var = "n")
                                   }
result_mvr_n_fac2_inter_df <- do.call("rbind",result_mvr_n_fac2_inter)
saveRDS(result_mvr_n_fac2_inter_df, "./linear_y_sim/mvr_n_fac2_inter.rds")

### AR1 rho
result_mvr_rho_lown250_ar1 <- foreach(i = rep(1:it,5), rho = rep(seq(0.1,0.9,0.2), each = it),
                                      .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba")) %dopar%{
                                        source('./RSource/mvr_julia.R')
                                        res <- sim_exp_mvr(n = 250, X.method = "AR1", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                           x.noise = 0, rho = rho, call.var = "rho")
                                      }
result_mvr_rho_lown250_ar1_df <- do.call("rbind",result_mvr_rho_lown250_ar1)
saveRDS(result_mvr_rho_lown250_ar1_df, "./linear_y_sim/mvr_rho_lown250_ar1.rds")


result_mvr_xnoise_fac2<- foreach(i = rep(1:it,4), xnoise = rep(c(0.1,0.5,1,2), each = it),
                                 .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                   source('./RSource/mvr_julia.R')
                                   res <- sim_exp_mvr(n = 500, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                      x.noise = xnoise, call.var = "x.noise")
                                 }
result_mvr_xnoise_fac2_df <- do.call("rbind",result_mvr_xnoise_fac2)
saveRDS(result_mvr_xnoise_fac2_df, "./linear_y_sim/mvr_xnoise.rds")

result_mvr_n_fac2_xnoise2 <- foreach(i = rep(1:it,11), N = rep(seq(150,250,10), each = it),
                                     .packages = c("mixOmics", "tidyverse", "knockoff", "CVglasso","rARPACK", "irlba", "JuliaCall")) %dopar%{
                                       source('./RSource/mvr_julia.R')
                                       res <- sim_exp_mvr(n = N, X.method = "fac2", y.noise = 0, p.in.cluster = 100, cluster = 5, seed.set = i,
                                                          x.noise = 2, rho = 0.5, call.var = "n")
                                     }
result_mvr_n_fac2_xnoise2_df <- do.call("rbind",result_mvr_n_fac2_xnoise2)
saveRDS(result_mvr_n_fac2_xnoise2_df, "./linear_y_sim/mvr_n_fac2_xnoise2.rds")
#result_mvr_n_fac2_df <- readRDS("./linear_y_sim/mvr_n_fac2.rds")

######## plotting ########
library(ggpubr)
result_mvr_s_fac2_df$s <- result_mvr_s_fac2_df$var
bench.s_fac2_fdr.plot <- ggplot(mapping = aes(x = s, y = fdp), data = result_mvr_s_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour = ko.method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "top")+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
bench.s_fac2_fdr.plot

bench.s_fac2_pow.plot <- ggplot(mapping = aes(x = s, y = tpp), data = result_mvr_s_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour= ko.method), fun.y="mean", geom="line") +
  theme(legend.position = "top")
bench.s_fac2_pow.plot

ggarrange(bench.s_fac2_fdr.plot, bench.s_fac2_pow.plot, ncol = 2, common.legend = T)


result_mvr_n_fac2_df$N <- result_mvr_n_fac2_df$var
bench.n_fac2_fdr.plot <- ggplot(mapping = aes(x = N, y = fdp), data = result_mvr_n_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour = ko.method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "top")+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
bench.n_fac2_fdr.plot

bench.n_fac2_pow.plot <- ggplot(mapping = aes(x = N, y = tpp), data = result_mvr_n_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour= ko.method), fun.y="mean", geom="line") +
  theme(legend.position = "top")
bench.n_fac2_pow.plot

ggarrange(bench.n_fac2_fdr.plot, bench.n_fac2_pow.plot, ncol = 2, common.legend = T)


result_mvr_noisey_fac2_df$y.noise <- result_mvr_noisey_fac2_df$var
bench.noisey_fac2_fdr.plot <- ggplot(mapping = aes(x = y.noise, y = fdp), data = result_mvr_noisey_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour = ko.method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "top")+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
bench.noisey_fac2_fdr.plot

bench.noisey_fac2_pow.plot <- ggplot(mapping = aes(x = y.noise, y = tpp), data = result_mvr_noisey_fac2_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour= ko.method), fun.y="mean", geom="line") +
  theme(legend.position = "top")
bench.noisey_fac2_pow.plot

ggarrange(bench.noisey_fac2_fdr.plot, bench.noisey_fac2_pow.plot, ncol = 2, common.legend = T)


###### caty ########
result_mvr_n_fac2_caty_df$N <- factor(result_mvr_n_fac2_caty_df$var, levels = as.character(seq(200, 1000, 100)))
bench.n_fac2_caty_fdr.plot <- ggplot(mapping = aes(x = N, y = fdp), data = result_mvr_n_fac2_caty_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour = ko.method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "top")+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
bench.n_fac2_caty_fdr.plot

bench.n_fac2_caty_pow.plot <- ggplot(mapping = aes(x = N, y = tpp), data = result_mvr_n_fac2_caty_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour= ko.method), fun.y="mean", geom="line") +
  theme(legend.position = "top")
bench.n_fac2_caty_pow.plot

ggarrange(bench.n_fac2_caty_fdr.plot, bench.n_fac2_caty_pow.plot, ncol = 2, common.legend = T)




######### quadratic #########
result_mvr_n_fac2_quad_df$N <- factor(result_mvr_n_fac2_quad_df$var, levels = as.character(seq(150, 500, 50)))
bench.n_fac2_quad_fdr.plot <- ggplot(mapping = aes(x = N, y = fdp), data = result_mvr_n_fac2_quad_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour = ko.method), fun.y="mean", geom="line") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "top")+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")
bench.n_fac2_quad_fdr.plot

bench.n_fac2_quad_pow.plot <- ggplot(mapping = aes(x = N, y = tpp), data = result_mvr_n_fac2_quad_df)+
  geom_boxplot(aes(fill = ko.method))+
  # scale_fill_manual(values = bench.col)+
  # scale_color_manual(values = bench.col)+
  stat_summary(aes(group= ko.method, colour= ko.method), fun.y="mean", geom="line") +
  theme(legend.position = "top")
bench.n_fac2_quad_pow.plot

ggarrange(bench.n_fac2_quad_fdr.plot, bench.n_fac2_quad_pow.plot, ncol = 2, common.legend = T)
