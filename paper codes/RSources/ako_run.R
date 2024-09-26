# AKO test
#setwd("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/linear_y_sim/ako")

source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/GraceKO_util.R")

AKO <- function(X, y, n_bootstraps = 25,  ko.method = "PLSKO.only", mvr = F,
                q = 0.05, offset = 0, w.method = "lasso.lcd", fdr.method = 'bhq', 
                gamma = 0.3, gamma_min = 0.1, seed = 1, ...){
  
  set.seed(seed)
  # n_bootstraps = 25
  pvals = matrix(0, ncol(X), n_bootstraps)
  selected <- list()
  #Multiple Knockoffs
  for (i in 1:n_bootstraps) {
    print(i)
    
    if(mvr){
      ko <- generate_ko_with_julia(X, Sigma = NULL, colMeans(X), method = "maxent", m = 1)$Xko
    }
    else ko <- knockoffX.sim(X, method = ko.method, ...)$ko
    
    S <- ko.filter(X = X, Xk = ko, y = y, q = q, method = w.method)
    
    pvals[,i] = empirical_pval(S$W, offset = offset)
    selected[i] <- S
  }
  
  aggregated_pval = apply(pvals, 1, quantile_aggregation, gamma=gamma, gamma_min=gamma_min, adaptive=FALSE)
  
  threshold = fdr_threshold(aggregated_pval, fdr=q, method = fdr.method,
                            reshaping_function=NULL)
  
  ako.s <- which(aggregated_pval <= threshold)
  
  return(list(t = threshold,
              w = aggregated_pval,
              s = selected,
              ako.s = ako.s
              ))
}

# test <- AKO(t(discv1.placenta.elevated$counts), discv1.placenta.elevated$samples$PE, w.method = "lasso.logistic", n_bootstraps = 50,
#             threshold.q = 0, r = 5, gamma = 0.5)
# which(test$w <= test$t)
# table(unlist(test$s))
# 
# test2 <- AKO(t(discv1.placenta.elevated$counts), discv1.placenta.elevated$samples$PE, w.method = "lasso.logistic", n_bootstraps = 50,
#              threshold.q = 0, r = 5, sparsity = 1, gamma = 0.3)
# saveRDS(test2, "plsko_50_b03.rds")
# which(test2$w <= test2$t)
# table(unlist(test2$s))
# 
# test3 <- AKO(t(discv1.placenta.elevated$counts), discv1.placenta.elevated$samples$PE, w.method = "lasso.logistic", n_bootstraps = 50,
#              threshold.q = 0, r = 5, sparsity = 0.8, gamma = 0.3)
# saveRDS(test3, "splsko_50_b03.rds")
# which(test3$w <= test3$t)
# table(unlist(test3$s))
# 
# 
# mvr.test <- AKO(t(discv1.placenta.elevated$counts), discv1.placenta.elevated$samples$PE, w.method = "lasso.logistic", n_bootstraps = 50,
#                 gamma = 0.3, mvr = T)
# mvr.test$ako.s
# saveRDS(mvr.test, "mvr_test_50_b03.rds")
# table(unlist(mvr.test$s))
# mvr.test$t

# default.ako <- AKO(t(discv1.placenta.elevated$counts), discv1.placenta.elevated$samples$PE, w.method = "lasso.logistic", n_bootstraps = 50, 
#              ko.method = "default")
# saveRDS(default.ako, "default_v19_ako.rds")
# default.ako$ako.s
