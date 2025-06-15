#' Generating multiple knockoff variables (MVR updated)
#'
#' This function samples multivariate Gaussian model-X knockoff variables multiple times for each original variable
#'
#' @param X A \code{n}-by-\code{p} matrix of original variables
#' @param mu A \code{p} vector of mean parameter of the Gaussian distribution of the original variables
#' @param Sigma A \code{p}-by-\code{p} covariance matrix of the original variables
#' @param omega A \code{p} vector indicating the weights for each variable. For now, we require each entry to be an integer greater than or equal to 2.
#' @param type Could be "entropy" (default), "sdp" or "equal", indicating the method that will be used to construct the knockoff variables.
#' @param diag_s A \code{p} vector of pre-computed (user-specified) covariance between the original and the knockoff variables. This will be computed according to method, if not supplied.
#'
#' @return a \code{n}-by-\code{p * max{omega - 1}} matrix containing the constructed knockoff variables. Although we construct the same number of knockoffs for all variables (which is the maximum cost), in subsequent steps, we only use the number of knockoffs based on the feature costs.
#'
#' @examples
#' library(cheapknockoff)
#' set.seed(123)
#' n <- 100
#' p <- 30
#' x <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(n)
#' omega <- c(2, 9, sample(seq(2, 9), size = 28, replace = TRUE))
#' # construct multiple knockoffs
#' X_k <- multiple_knockoff_Gaussian(X = x, mu = rep(0, p), Sigma = diag(1, p), omega = omega)
#' @import CVXR
#' @export
multiple_knockoff_Gaussian <- function(X, num_knockoff, mu = colMeans(X), Sigma = cov(X), type = c("entropy", "sdp", "equi", "mvr"), diag_s = NULL){
  # generate multiple knockoff variables of the orignal X, which
  # follows a Gaussian(mu, Sigma)
  # X_j is constructed omega_j - 1 times for each j

  start.time <- Sys.time()
  type <- match.arg(type)

  n <- nrow(X)
  p <- ncol(X)
  
  # first compute the diag_s, which is the diagonal matrix
  # that captures the difference in covariance between the orignal and the knockoffs
  if (is.null(diag_s)){
    diag_s = diag(as.numeric(switch(match.arg(type),
                                    entropy = solve_entropy(Sigma, num_knockoff),
                                    sdp = solve_sdp(Sigma, num_knockoff),
                                    equi = solve_equi(Sigma, num_knockoff),
                                    mvr = solve_mvr(Sigma, num_knockoff))))
  }
  
  # after the diagonal matrix is computed
  # we can generate the Gaussian multiple knockoffs by
  # sampling from the conditional distribution given the original variables
  # which is still a Gaussian
  
  # intermediate variable: Sigma^{-1} D
  Siginv_diag_s <- solve(Sigma, diag_s)
  
  # mean vector for the k-th knockoff (1 <= k <= num_knockoff)
  mu_k <- matrix(rep(as.numeric(crossprod(Siginv_diag_s, mu)), n), nrow = n, byrow = TRUE) + X %*% t(diag(rep(1, p)) - Siginv_diag_s)
  # mean vector is num_knockoff copies of mu_k
  mu_full <- matrix(rep(1, num_knockoff), nrow = 1) %x% mu_k
  
  # compute the conditional convariance matrix C
  CmD <- diag_s - diag_s %*% Siginv_diag_s
  # here for simplicity we use Kronecker product
  C <- matrix(1, nrow = num_knockoff, ncol = num_knockoff) %x% CmD + diag(rep(1, num_knockoff)) %x% diag_s
  
  stopifnot(isSymmetric(C))
  
  # the Cholesky decomposition of C
  # recall C_chol^T C_chol = C
  eps <- 1e-4
  min_lam <- min(eigen(C, only.values = TRUE)$values)
  if(min_lam <= eps){
    print('Generated diagonal is invalid. Adding eps to the diagonal of the knockoff joint covariance matrix.')
    C <- C + diag(abs(min_lam) + eps, dim(C))
  }
  
  # Finally generate the desired knockoff variables
  x <- matrix(rnorm(n * p * num_knockoff), nrow = n, ncol = p * num_knockoff)
  x <- mu_full + x %*% chol(C)
  
  return(list(X_k = x, run.time = difftime(Sys.time(), start.time, units = "secs")))
}


# call CVXR to get an estimated of diagonal matrix using Entropy Maximization
solve_entropy <- function(Sigma, num_knockoff, eps = 1e-8){
  Sig_scaled <- (num_knockoff + 1) / num_knockoff * Sigma
  # ambient dimension p
  p <- ncol(Sigma)
  # construct optimization variable D = (s_j)
  D <- Variable(p)
  # objective definition
  obj <- Maximize(log_det(Sig_scaled - diag(D)) + num_knockoff * sum_entries(log(D)))
  # constraint definition
  # CVXR Version 1.0 Note: Positive semi-definite variables are now designated using PSD = TRUE rather than the Semidef function!
  if (packageVersion("CVXR") > "0.99-7") {
    S  <- Variable(p, p, PSD = TRUE)
    constr <- list(Sig_scaled - diag(D) - diag(eps, p) == S, D >= 0)
  } else {
    constr <- list(Sig_scaled - diag(D) - diag(eps, p) == Semidef(p), D >= 0)
  }
  # equivalently
  # constr <- list(lambda_min(Sig_scaled - diag(D)) > 0,
  #                D >= 0)
  
  prob <- Problem(objective = obj, constraints = constr)
  result <- psolve(prob)
  # if(result$status != "optimal")
  #   cat("CVX does not fully solve the optimization problem, the status is:", result$status, fill = TRUE)
  return(result$getValue(D))
}

# call CVXR to get an estimated of diagonal matrix using SDP
# solve_sdp <- function(Sigma, num_knockoff){
#   eps <- 1e-4
#   Sig_scaled <- (num_knockoff + 1) / num_knockoff * Sigma
#   # ambient dimension p
#   p <- ncol(Sigma)
#   # construct optimization variable D = (s_j)
#   D <- Variable(p)
#   # objective definition
#   obj <- Minimize(sum_entries(abs(rep(1, p) - D)))
#   # constraint definition
#   #constr <- list(Sig_scaled - diag(D) - diag(eps, p) == Semidef(p),
#   #               D >= 0)
#   # or equivalently
#   constr <- list(lambda_min(Sig_scaled - diag(D)) > 0,
#                  D >= 0)
#   
#   prob <- Problem(objective = obj, constraints = constr)
#   result <- psolve(prob)
#   # if(result$status != "optimal")
#   #   cat("CVX does not fully solve the optimization problem, the status is:", result$status, fill = TRUE)
#   return(result$getValue(D))
# }

# call CVXR to get an estimated Equi-correlated Knockoffs
solve_equi <- function(Sigma, num_knockoff){
  sig_min <- min(eigen(Sigma, only.values = TRUE)$value) * (num_knockoff + 1) / num_knockoff
  # ambient dimension p
  p <- ncol(Sigma)
  # result is just the diagonal matrix with entries having same values of sig_min
  return(rep(sig_min, p))
}


solve_mvr <- function(Sigma, num_knockoff, max_iter = 10000){
  # scale Sigma
  Sig_scaled <- (num_knockoff + 1) / num_knockoff * Sigma
  
  p <- dim(Sigma)[1]
  I <- diag(p)
  # Inital value
  S <- diag((min(eigen(Sigma, symmetric = T)$value)), p)
  D <- Sig_scaled - S
  L <- t(chol(D))
  
  for (iter in 1:max_iter){
    S.old <- S
    for (j in 1:p) {
      v_d <- forwardsolve(l = L, x = I[,j])
      v_n <- backsolve(r = t(L), x = v_d)
      c_d <- sum(v_d^2)
      c_n <- -sum(v_n^2)
      
      delta <- optimise(function(x){1/(S[j,j]+x)-x*c_n/(1-x*c_d)},
                        lower = -S[j,j],
                        upper = 1/c_d,
                        maximum = F)$minimum
      
      S[j,j] <- S[j,j] + delta
      L <- t(chol(Sig_scaled - S))
    }
    if(max(abs(S - S.old)) < 1e-13){
      break()
    }
  }
  S <- matrix(diag(S))
  return(S)
}

solve_sdp <- function(Sigma, num_knockoff, gaptol=1e-6, maxit=1000, verbose=FALSE) {
  # Check that covariance matrix is symmetric
  stopifnot(isSymmetric(Sigma))
  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]
  
  scale = (num_knockoff + 1) / num_knockoff
  
  # Check that the input matrix is positive-definite
  if (!is_posdef(G)) {
    warning('The covariance matrix is not positive-definite: knockoffs may not have power.', immediate.=T)
  }
  
  # Convert problem for SCS
  
  # Linear constraints
  Cl1 = rep(0,p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1,p)
  Al2 = Matrix::Diagonal(p)
  
  # Positive-definite cone
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x=d_As)
  As = As[which(Matrix::rowSums(As) > 0),] 
  Cs = c(scale*G)
  
  # Assemble constraints and cones
  A = cbind(Al1,Al2,As)
  C = matrix(c(Cl1,Cl2,Cs),1)
  K=NULL
  K$s=p
  K$l=2*p
  
  # Objective
  b = rep(1,p)
  
  # Solve SDP with Rdsdp
  OPTIONS=NULL
  OPTIONS$gaptol=gaptol
  OPTIONS$maxit=maxit
  OPTIONS$logsummary=0
  OPTIONS$outputstats=0
  OPTIONS$print=0
  if(verbose) cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A,b,C,K,OPTIONS)
  if(verbose) cat("done. \n")
  
  # Check whether the solution is feasible
  if( ! identical(sol$STATS$stype,"PDFeasible")) {
    warning('The SDP solver returned a non-feasible solution. Knockoffs may lose power.')
  }
  
  # Clip solution to correct numerical errors (domain)
  s = sol$y
  s[s<0]=0
  s[s>1]=1
  
  # Compensate for numerical errors (feasibility)
  if(verbose) cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-8
  while ((psd==0) & (s_eps<=0.1)) {
    if (is_posdef(scale*G-diag(s*(1-s_eps),length(s)),tol=1e-9)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)
  s[s<0]=0
  if(verbose) cat("done. \n")
  
  # Verify that the solution is correct
  if (all(s==0)) {
    warning('In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.',immediate.=T)
  }
  
  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}
