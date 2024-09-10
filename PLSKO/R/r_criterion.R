#' Calculate the Empirical Number of Components for PLS Regression
#'
#' This function calculates the empirical number of components used for PLS regression
#' by applying the \eqn{PC_p1} criterion from Bai and Ng (2002). It determines the optimal number of factors based on the minimum value of the criterion. The method and function is adapted from Y Fan et al. (2019).
#'
#' @param X A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param rmax An integer specifying the maximum number of factors to consider. Default is 10.
#'
#' @return An integer representing the optimal number of components based on the \eqn{PC_p1} criterion.
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100*10), 100, 10)
#' optimal_r <- r_criterion(X)
#'
#' @references Fan Y, Lv J, Sharifvaghefi M et al. IPAD: Stable Interpretable Forecasting with Knockoffs Inference. Journal of the American Statistical Association 2020;115:1822–34.
#' @references Bai J, Ng S. Determining the Number of Factors in Approximate Factor Models. Econometrica 2002;70:191–221.


#' @export
#'
r_criterion <- function(X, rmax = 10){
  # Number of observations (n) and variables (p)
  n = nrow(X)
  p = ncol(X)

  # Initialize a matrix to store the PC_p1 criterion values for each possible number of components
  PC = matrix(NA,rmax+1,1)

  SX = scale(X)
  mXX = SX %*%t(SX)

  # Calculate the criteria through the number of components from rmax to 0
  for(k in rmax:0){
    if (k == 0){
      PC[k+1,] = sum(SX^2/(n*p))
    } else
      eig <- RSpectra::eigs_sym(mXX, k, which = "LM", sigma = NULL, lower = TRUE)
    meigvec <- eig$vectors
    mF <- sqrt(n) * meigvec       # estimated factors
    Lam <- (t(mF) %*% SX)/n
    if(k==rmax){
      sigma2 = sum((SX-mF %*% Lam)^2)/(n*p)
    }
    PC[k+1,] = sum((SX-mF %*% Lam)^2)/(n*p) + k*sigma2*((n+p)/(n*p))*log((n*p)/(n+p))
  }
  r = which(PC == min(PC))[1]-1
}
