
# PART A: KNOCKOFFS WITH JULIA (see also https://github.com/biona001/Knockoffs.jl)
        
# PART B: KNOCKOFFS FOR E-MLKF (see also https://katsevich-lab.github.io/publications/)

# PART C: CALCULATE KNOCKOFF E-VALUES FOR UKB + MISC KNOCKOFF FUNCTIONS


#### PART A: KNOCKOFFS WITH JULIA ########

# This knockoff construction is based on "Variable Selection with Knockoffs" by Benjamin Chu 
# See https://github.com/biona001/Knockoffs.jl

generate_ko_with_julia <- function(X, Sigma, mu, method, m) {
  # pass variables from R to Julia
  julia_assign("X", X)
  julia_assign("Sigma", Sigma)
  julia_assign("mu", mu)
  julia_assign("method", method)
  julia_assign("m", m)
  # with Julia, solve group knockoffs problem
  julia_command("result = modelX_gaussian_knockoffs(X, Symbol(method), mu, Sigma, m=Int(m), verbose=true)")
  
  # pull variables from Julia back to R
  Xko <- julia_eval("result.X̃ ") # the knockoffs
  # return
  result <- list("Xko"=Xko)
  return (result)
}

run_group_ko_lasso_with_julia <- function(y, X, Sigma, mu, method, m, groups,alpha) {
  # pass variables from R to Julia
  julia_assign("y", y)
  julia_assign("X", X)
  julia_assign("Sigma", Sigma)
  julia_assign("mu", mu)
  julia_assign("method", method)
  julia_assign("m", m)
  julia_assign("groups", groups)
  julia_assign("fdrs", alpha)
  
  # with Julia, generate group knockoffs then solve lasso
  julia_command("result = fit_lasso(vec(y), X, mu, Symmetric(Sigma), groups=Int.(groups), method=Symbol(method), m=Int(m), fdrs=[fdrs])")
  
  # pull variables from Julia back to R
  selected <- julia_eval("result.selected") # selected variables
  Xko <- julia_eval("result.ko.X̃ ") # the knockoffs
  S <- julia_eval("result.ko.S") # the S matrix
  obj <- julia_eval("result.ko.obj") # final objective value
  W <- julia_eval("result.W")
  taus <- julia_eval("result.τs")
  fdr_target <- julia_eval("result.fdr_target")
  
  # return
  result <- list("selected"=selected, "Xko"=Xko, "S"=S, "obj"=obj, "W" = W, "taus" = taus, "fdr_target" = fdr_target)
  return (result)
}

run_ko_lasso_with_julia <- function(y, X, Sigma, mu, method, m,alpha) {
  # pass variables from R to Julia
  julia_assign("y", y)
  julia_assign("X", X)
  julia_assign("Sigma", Sigma)
  julia_assign("mu", mu)
  julia_assign("method", method)
  julia_assign("m", m)
  julia_assign("fdrs", alpha)
  
  # with Julia, generate group knockoffs then solve lasso
  julia_command("result = fit_lasso(vec(y), X, mu,  Symmetric(Sigma), method=Symbol(method), m=Int(m), fdrs=[fdrs])")
  
  # pull variables from Julia back to R
  selected <- julia_eval("result.selected") # selected variables
  Xko <- julia_eval("result.ko.X̃ ") # the knockoffs
  W <- julia_eval("result.W")
  taus <- julia_eval("result.τs")
  fdr_target <- julia_eval("result.fdr_target")
  
  
  # return
  result <- list("selected"=selected, "Xko"=Xko, "W" = W, "taus" = taus, "fdr_target" = fdr_target)
  return (result)
}




############ PART B: E-MLKF KNOCKOFFS####################


# See Katsevich, Sabatti (2019) Multilayer knockoff filter: Controlled variable selection at multiple resolutions

# Function downloaded from: https://katsevich-lab.github.io/publications/

# compute knockoff variables at each layer
get_knockoff_variables = function(X, groups, knockoff_type){
  Sys.sleep(1)
  set.seed(2022)
  Sys.sleep(1)
  switch(knockoff_type,
         fixed_equi = {
           X.knockoffs = list()
           for(m in 1:M){
             cat(sprintf("Constructing knockoff variables for layer %d...\n", m))
             X.knockoffs[[m]] = construct_knockoffs_equi(X, groups[,m])
           }
         }       
  )
  return(X.knockoffs)  
}


# Function downloaded from: https://katsevich-lab.github.io/publications/

# construct equicorrelated fixed-design knockoffs
construct_knockoffs_equi = function(X, groups){
  
  Sys.sleep(1)
  set.seed(2022)
  Sys.sleep(1)
  
  # problem dimensions
  N = nrow(X)
  n = ncol(X)
  stopifnot(N >= 2*n) # make sure problem is low-dimensional
  
  # make sure columns of X are normalized to have norm 1
  Sigma = crossprod(X)
  stopifnot(max(abs(diag(Sigma)-1)) <  1e-10)
  
  #print(summary(rowMeans(Sigma)))
  
  # equicorrelated construction
  D = matrix(0, n, n)
  S = matrix(0, n, n)
  G = max(groups)
  for(g in 1:G){
    group_idx = which(groups == g)
    group_size = length(group_idx)
    Sigma_grp = Sigma[group_idx, group_idx]
    S[group_idx, group_idx] = Sigma_grp
    eig = eigen(Sigma_grp)
    D[group_idx, group_idx] = eig$vectors%*%diag(1/sqrt(eig$values), 
                                                 group_size, group_size)%*%t(eig$vectors)
  }
  
  gamma = min(1, 2*min(eigen(D%*%Sigma%*%D, symmetric = TRUE)$values)) - 1e-10
  S = gamma*S

  set.seed(2022)
  U_tilde = eigen(crossprod(t(X)))$vectors[,(n+1):(2*n)]
  C = chol(2*S - S%*%solve(Sigma,S))
  X.knockoff = X%*%(diag(n) - solve(Sigma, S)) + U_tilde%*%C
 
  return(X.knockoff)
}

# Function is downloaded from: https://katsevich-lab.github.io/publications/; Katsevich, Sabatti (2019) Multilayer knockoff filter: Controlled variable selection at multiple resolutions
# INPUT: 
  # X
  # X.knockoff
  # y 
  # groups: data frame of dimensions p x M (n number of individual hypotheses and M the number of resolutions) 
# OUTPUT: 
  # W: importance statistics

# get group lasso signed max statistics 
get_LSM_stats = function(X, X.knockoff, y, groups){
  
  set.seed(2022)
  
  # get problem dimensions
  n = ncol(X)
  G = max(groups)
  
  # make sure that groups are in consecutive order for gglasso
  ord = order(groups)
  X = X[,ord]
  X.knockoff = X.knockoff[,ord]
  groups = groups[ord]
  
  # augment the problem  
  groups_augmented = c(groups, groups + G)
  X_augmented = cbind(X, X.knockoff) 
  
  # compute lambda_max  
  inner_prods = t(X_augmented)%*%y
  inner_prods_sq = inner_prods*inner_prods
  inner_product_norms = sqrt(sapply(1:(2*G), function(g)(mean(inner_prods_sq[groups_augmented == g]))))
  lambda_max = max(inner_product_norms)/nrow(X)
  
  # define lambda grid  
  num_lambda = 1000
  lambda_values = seq(lambda_max, 0, length = num_lambda)
  
  # run group lasso, using glmnet if groups all have size one because it's faster
  if(n > G){
    output = gglasso(X_augmented, y, group = groups_augmented,
                     lambda = lambda_values, intercept=FALSE)
  } else if(n == G){
    output = glmnet(X_augmented, y, lambda = lambda_values, standardize = FALSE, intercept = FALSE)
  }
  
  # compute first entry times  
  min_idx = sapply(1:(2*n), function(row)(min(c(num_lambda, which(output$beta[row,] != 0)))))
  first_entries = sapply(1:(2*G), function(g)(min(min_idx[groups_augmented == g])))
  max_lambdas = lambda_values[first_entries]
  
  # compute knockoff statistics
  W = sapply(1:G,
             function(g)(max(max_lambdas[g],
                             max_lambdas[g+G])*sign(max_lambdas[g] - max_lambdas[g+G])))
  
  return(W)
}



# Function is downloaded from: https://katsevich-lab.github.io/publications/; Katsevich, Sabatti (2019) Multilayer knockoff filter: Controlled variable selection at multiple resolutions
# INPUT: 
  # X
  # X.knockoff
  # y 
  # groups: data frame of dimensions p x M (n number of individual hypotheses and M the number of resolutions) 

# OUTPUT: 
  # W: importance statistics

# compute knockoff statistics at each layer
get_knockoff_stats = function(X, X.knockoffs, y, groups, statistic_type){
  #source("aux_knockoffs.R")
  
  set.seed(2022)
  
  # problem dimensions
  n = nrow(groups)
  M = ncol(groups)
  G = apply(groups,2,max) 
  
  # initialize
  W = list()
  
  # compute W layer by layer
  for(m in 1:M){
    cat(sprintf("Computing knockoff statistics at layer %d...\n", m))
    # get knockoffs at this layer
    X.knockoff = X.knockoffs[[m]]
    # switch based on statistic type
    switch(statistic_type,
           # group lasso signed max (Fisher-like)
           group_LSM = {
             reg_groups = groups[,m]
             W[[m]] = get_LSM_stats(X, X.knockoff, y, reg_groups)
           },
           
           # individual lasso signed max (Simes-like)
           ind_LSM = {
             reg_groups = 1:n
             W_ind = get_LSM_stats(X, X.knockoff, y, reg_groups)
             FUN = function(g){
               if(max(W_ind[groups[,m] == g]) + min(W_ind[groups[,m] == g]) == 0){
                 return(0)
               } else{
                 return(W_ind[groups[,m] == g][which.max(abs(W_ind[groups[,m] == g]))])
               }
             }
             W[[m]] = sapply(1:G[m], FUN)
           }
    )
  }
  return(W)
}


#### PART C: KNOCKOFF E-VALUES + MISC FUNCTIONS ########

# function to calculate partial conjunction e-values, uses the function partial_conjunction_averaging below
# INPUT: 
# u: partial conjunction u to be used 
# W: importance statistics
# Q: number of outcomes 
# M: number of resolutions 
# gamma: gamma to be used in e-value definition 
# group_size_with_individual: vector specifying the group sizes in each layer, including individual level with group size = 1
# offset: offset to be used for knockoff threshold
# OUTPUT: list containing
# evals_global_per_m: list of e-values for each level of resolution 
# fracs_global_per_m: list of e-values without multiplicative factor of group size ("fraction") for each level of resolution
calc_evalues_partial <- function(u = 1, W, Q, M, p, gamma, alpha, group_size_with_individual, offset = 1) {
  
  evals_global_per_m <- list()
  fracs_global_per_m <- list()
  global_rejected_per_layer <- list() 
  taus_evals <- list()
  
  all_evals <- list()
  
  for(m in 1:M) {
    
    ngroups_in_m <- p / group_size_with_individual[m]
    evals_within_m <- matrix(NA, nrow = ngroups_in_m , ncol = Q)
    fracs_within_m <- matrix(NA, nrow = ngroups_in_m , ncol = Q)
    taus_evals_m <- list()
    
    
    for(q in 1:Q) {
      
      # no replacement of gamma
      taus_evals_m[[q]] = knockoff.threshold(W[[q]][[m]], fdr=gamma, offset=offset)
      
      
      evals_within_m[, q] <- ngroups_in_m * ((W[[q]][[m]] >=  taus_evals_m[[q]] ) / (1 + sum(W[[q]][[m]] <= -taus_evals_m[[q]] )))
      fracs_within_m[, q] <- ((W[[q]][[m]] >=  taus_evals_m[[q]] ) / (1 + sum(W[[q]][[m]] <= -taus_evals_m[[q]] )))
      
    }
    
    evals_global_per_m[[m]] <- partial_conjunction_averaging(evals_within_m, Q = Q, u = u)$avg
    fracs_global_per_m[[m]] <- partial_conjunction_averaging(fracs_within_m, Q = Q, u = u)$avg
    
  }
  
  
  
  return(list(evals = evals_global_per_m, fracs = fracs_global_per_m ))
  
}


# function to calculate partial conjunction e-values (ie. average the Q - u + 1 smallest)
# INPUT: 
# df: data frame containing the e-values 
# Q: number of outcomes (or number of individual hypotheses considered in the partial conjunction null)
# u: partial conjunction u

# OUTPUT: 
# partial_conjunction_avg: partial conjunction e-values
partial_conjunction_averaging <- function(df, Q, u) {
  df <- data.frame(df)
  colnames(df) <- paste0("Y", seq(1:Q))
  
  df %>% 
    rownames_to_column(var = "feature") %>%
    pivot_longer(cols = -feature) %>%
    mutate(feature = as.integer(feature)) %>% 
    group_by(feature) %>%
    arrange(value, .by_group = TRUE) %>%
    slice(1:(Q - u + 1)) %>% 
    summarise(avg = mean(value)) -> partial_conjunction_avg
  
  return(partial_conjunction_avg)
}


# code modified from: https://github.com/cran/knockoff/blob/master/R/knockoff.R; Candes et al., "Panning for gold: model-X knockoffs for high-dimensional controlled variable selection", J. R. Statist. Soc. B (2018) 80, 3, pp. 551-577.

# function to find the minimum possible "alpha" for which there would be at least one rejection
# INPUT: 
# W: importance statistics 
# offset: offset to be used in knockoff 
# OUTPUT: a list containing: 
# the threshold corresponding to the minimum alpha 
# the minimum possible alpha
min.alpha.knockoff.threshold <- function(W, offset=1) {
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  min_alpha <- min(ratio)
  threshold <- ts[which.min(ratio)]
  return(list(threshold, min_alpha))
}



