AKO <- function(X, y, n_bootstraps = 25,  ko.method = "PLSKO", mvr = F,
                q = 0.05, offset = 1, w.method = "lasso.lcd",
                gamma = 0.3, seed = 1, ...){
  
  set.seed(seed)
  # n_bootstraps = 25
  pvals = matrix(0, ncol(X), n_bootstraps)
  selected <- list()
  #Multiple Knockoffs
  for (i in 1:n_bootstraps) {
    # if(mvr){
    #   ko <- generate_ko_with_julia(X, Sigma = NULL, colMeans(X), method = "maxent", m = 1)$Xko
    # }
    # else ko <- knockoffX.sim(X, method = ko.method, ...)$ko
    ko = plsko(X, ...)$X_k
    S <- ko.filter(X = X, Xk = ko, y = y, q = q, method = w.method)

    pvals[,i] = empirical_pval(S$W, offset = offset)
    selected[i] <- S
  }
  
  aggregated_pval = apply(pvals, 1, quantile_aggregation, gamma=gamma)
  
  threshold = bhq_threshold(aggregated_pval, fdr=q)
  
  ako.s <- which(aggregated_pval <= threshold)
  
  return(list(s = selected,
              ako.s = ako.s
  ))
}

# AKO with generated knockoff sets as a list
AKO_with.ko <- function(X, y, Xko.list,
                    q = 0.05,  w.method = "lasso.lcd", offset = 1,
                    gamma = 0.3, seed = 1){
  
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
  
  aggregated_pval = apply(pvals, 1, quantile_aggregation, gamma=gamma)
  
  threshold = bhq_threshold(aggregated_pval, fdr=q)
  
  ako.s <- which(aggregated_pval <= threshold)
  
  return(list(s = selected,
              ako.s = ako.s
  ))
}

empirical_pval = function(test_score, offset = 1){
  pvals = c()
  n_features = length(test_score)
  if (offset !=0 && offset!=1){
    return("'offset' must be either 0 or 1")
  }
  else{
    test_score_inv = -test_score
    for (i in 1:n_features){
      if (test_score[i] <= 0){
        pvals = c(pvals, 1)
      }
      else{
        
        pvals = c(pvals,(offset+sum(test_score_inv[i] >= test_score))/n_features)
      }
    }
  }
  return (pvals)
}

quantile_aggregation = function(pvals, gamma=0.3){
  converted_score = (1 / gamma) *  quantile(pvals, gamma)
  return (min(1, converted_score))
}

bhq_threshold = function(pvals, fdr=0.1){
  n_features = length(pvals)
  pvals_sorted = sort(pvals)
  selected_index = 2 * n_features
  for (i in seq(n_features, 1, -1)){
    if (pvals_sorted[i] <= (fdr * i / n_features)){
      selected_index = i
      break
    }
  }
  if (selected_index <= n_features){
    return (pvals_sorted[selected_index])
  }
  else{
    return ('-1.0')
  }
}
