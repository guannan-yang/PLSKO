# AKO with generated knockoff variable set as list
 
AKO.wko <- function(X, y, Xko.list,
                q = 0.05,  w.method = "lasso.lcd", offset = 0, fdr.method = 'bhq', 
                gamma = 0.3, gamma_min = 0.1, seed = 1, ...){
  
  set.seed(seed)
  
  n_bootstraps = length(Xko.list)
  pvals = matrix(0, ncol(X), n_bootstraps)
  selected <- list()
  #Multiple Knockoffs
  for (i in 1:n_bootstraps) {
    print(i)
    
    ko <- Xko.list[[i]]
    
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

### top variables 200 ====
# most.variable <- discv1$counts[rank(apply(discv1$counts, 1, sd))>(nrow(discv1$counts)-200),]
# 
# discv1.filter <- discv1
# discv1.filter$counts <- most.variable
# 
# discv1.filter$samples$PE <- ifelse(discv1.filter$samples$disease=="control", yes = "control", no="PE")
# discv1.filter$samples$PE <- as.factor(discv1.filter$samples$PE)
# 
# discv1_plsko_200 <- readRDS("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/RealData/Stanford Nature/analysis/discv1_plsko_200.rds")
# discv1_plsko_200_ko <- lapply(discv1_plsko_200, function(x){
#   x$ko$ko
# })
# 
# test.200 <- AKO.wko(t(discv1.filter$counts), y = discv1.filter$samples$PE, Xko.list = discv1_plsko_200_ko, 
#                     w.method = "lasso.logistic")
# saveRDS(test.200, "ako_200.rds")
# table(unlist(test.200$s))


### top 200 semi-sim ====
# set.seed(1)
# y.200 <- DGPy.AR1(t(discv1.filter$counts), s = 25, A = c(3,5))
# test.200 <- AKO.wko(t(discv1.filter$counts), y = y.200$y, Xko.list = discv1_plsko_200_ko, 
#                     w.method = "lasso.lcd")
# fdr(test.200$ako.s, y.200$Beta)
# pow(test.200$ako.s, y.200$Beta)

# > fdr(test.200$ako.s, y.200$Beta)
# [1] 0.1764706
# > pow(test.200$ako.s, y.200$Beta)
# [1] 0.56

# set.seed(1)
# y.200 <- DGPy.AR1(t(discv1.filter$counts), s = 25, A = c(3,5))
# test.200.plus <- AKO.wko(t(discv1.filter$counts), y = y.200$y, Xko.list = discv1_plsko_200_ko, 
#                     w.method = "lasso.lcd", offset = 1)
# fdr(test.200.plus$ako.s, y.200$Beta)
# pow(test.200.plus$ako.s, y.200$Beta)