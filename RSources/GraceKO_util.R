## #' Function of intermediary products derived from the conversion of feature statistics (single knockoffs) to aggregated feature statistics.
#'
#' @param test_score a vector of the feature statistics obatianed from single knockoff procedure
#' @param offset a modified value of the conversion procedure
#'
#' @return a vector of the intermediary products

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

quantile_aggregation = function(pvals, gamma=0.5, gamma_min=0.00001, adaptive=FALSE){
  if (adaptive == TRUE){
    return (adaptive_quantile_aggregation(pvals, gamma_min))
  }
  
  else{
    return (fixed_quantile_aggregation(pvals, gamma))
  }
  
}

#' Function to apply quantile aggregation function based on Meinshausen et al (2008)
#'
#' @param pvals a n_bootstrap by n_covariates matrix of the feature statistics, AGCs
#' @param gamma the pre-specified percentile value used for aggregation, range from 0 to 1

fixed_quantile_aggregation = function(pvals, gamma = 0.3){
  
  
  converted_score = (1 / gamma) *  quantile(pvals, gamma)
  
  return (min(1, converted_score))
}



#' Function to apply adaptive version of the quantile aggregation method, Meinshausen et al. (2008)
#'
#' @param pvals a n_bootstrap by n_covariates matrix of the feature statistics, AGCs
#' @param gamma_min the pre-specified percentile value used for aggregation, range from 0 to 1

adaptive_quantile_aggregation = function(pvals, gamma_min=0.05){
  
  gammas = seq(gamma_min, 1.05, 0.05)
  list_Q = matrix(0,nrow = length(gammas), ncol = length(gammas))
  for (gamma in gammas) {
    list_Q = c(list_Q, fixed_quantile_aggregation(pvals, gamma))
    
  }
  
  return (min(1, (1 - log(gamma_min)) * min(list_Q)))
  
}

#' Function to compute the corresponding threshold (BY/BH) for controlling FDR
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param fdr the pre-specified false discovery rate level
#' @param method the test procedure applied to compute the threshold. The defult method to compute threshold is the the Standard Benjamini-Hochberg procedure.
#' @param reshaping_function the default value for reshaping function can be referred by Benjamini & Yekutieli (2001) otherwise, the reshaping function can be referred by Ramdas et al (2017)
#'
#' @return The threshold value corresponding to the test procedure

fdr_threshold = function(pvals, fdr=0.1, method='bhq', reshaping_function=NULL){
  if (method == 'bhq'){
    # pvals_bhq = as.vector(unlist(pvals))
    return (bhq_threshold(pvals, fdr=fdr))
  }
  else{
    if(method == 'bhy'){
      # pvals_bhy = as.vector(unlist(pvals))
      return( bhy_threshold(
        pvals, fdr=fdr, reshaping_function=reshaping_function))
    }
    else{
      return('{} is not support FDR control method')
    }
  }
}


#' Function for applying the Standard Benjamini-Hochberg for controlling False discovery rate
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param fdr the pre-specified false discovery rate level
#'
#' @return The threshold value corresponding to the Standard Benjamini-Hochberg procedure

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

#' Function to apply the Benjamini-Hochberg-Yekutieli procedure for controlling FDR
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param reshaping_function the default value for reshaping function can be referred by Benjamini & Yekutieli (2001) otherwise, the reshaping function can be referred by Ramdas et al (2017)
#' @param fdr the pre-specified false discovery rate level
#' @return The threshold value corresponding tthe Benjamini-Hochberg-Yekutieli procedure


bhy_threshold = function(pvals,
                         reshaping_function = NULL, fdr = 0.1)
{
  
  n_features = length(pvals)
  p_vals_sorted = sort(pvals)
  selected_index = 2 * n_features
  
  
  if (is.null(reshaping_function)==TRUE)
  {
    temp = seq(n_features,1)
    sum_inverse = 0
    for (i in temp) {
      sum_inverse = sum_inverse + 1 / i
    }
    
    return (bhq_threshold(pvals, fdr / sum_inverse))
  }
  else{
    for (i in seq(n_features - 1, 0, -1)){
      if (p_vals_sorted[i] <= fdr * reshaping_function(i) / n_features){
        selected_index = i
        break
      }
      
    }
    
    if (selected_index <= n_features){
      return (p_vals_sorted[selected_index])
    }
    else{
      return ('-1.0')
    }
    
  }
  
}
