#' Pipeline for the whole knockoff framework with PLSKO generated knockoff variables
#'
#' @import progress
#' @param X A numeric matrix or data frame of predictors.
#' @param y A numeric vector or factor of responses.
#' @param q A numeric value specifying the false discovery rate (FDR) level. Default is `0.05`.
#' @param method A string specifying the method to compute test statistics. Options are `"lasso.lcd"`, `"lasso.logistic"`, `"lasso.max.lambda"`, `"RF"`. Default is `"lasso.lcd"`.
#' @param offset A numeric or string value to adjust the knockoff threshold. Default is `0` (control modifed FDR). Other options include `1` (yields a slightly more conservative procedure ("knockoffs+") that controls the FDR according to the usual definition) and `"both"` (returns results from both "knockoffs" and "knockoffs+").
#' @param nb.list Optional. A list of length \eqn{p} or adjacency matrix of \eqn{p \times p} that defines the neighbourship of variables. A list of length \eqn{p} should include the neighbours' index of each variable from \eqn{X_1} to \eqn{X_p} in order; The \eqn{i^{th}} element in the list includes the indices of the neighbour variables of \eqn{X_i}, or \code{NULL} when no neighbours. A adjacency matrix should be symmetric with only binary element and  where \eqn{M_{ij} = 1} when \eqn{X_i} and \eqn{X_j} is defined as neighbours; otherwise \eqn{M_{ij} = 0} when not neighbour or on diagonal (i.e. \eqn{i = j}). If not provided or NULL, the neighborhoods are determined based on correlations.
#' @param threshold.abs Optional. A value between \eqn{0} and \eqn{1}. A numeric value specifying an absolute correlation threshold to define neighborhoods. If not provided, the function uses \code{threshold.q}.
#' @param threshold.q Optional. A numeric value between 0 and 1 indicating the quantile of the correlation values to use as a threshold.
#' @param ncomp Optional. An integer specifying the number of components to use in the PLS regression. Default is \code{NULL}, the \code{ncomp} is determined empirically by \eqn{PC_p1} criterion.
#' @param sparsity Optional. A numeric value between 0 and 1 specifying the sparsity level in the PLS regression. Default is 1 (no sparsity).
#' @param rmax An integer specifying the maximum number of factors to consider when \code{ncomp} is not defined. Default is 5.
#'
#' @param seed An integer seed for reproducibility. Default is 1.
#'
#' @return An object of class "knockoff.result". This object is a list containing at least the following components:
#' \describe{
#'   \item{`call`}{The matched call.}
#'   \item{`statistic`}{The computed test statistics `W` of length \eqn{p}.}
#'   \item{`threshold`}{The knockoff threshold used to determine discoveries.}
#'   \item{`selected`}{The named vector of selected variables.}
#' }
#' @export
#'
#' @examples
#' # Example 1: continuous response
#' set.seed(1)
#' X <- matrix(rnorm(100*10), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' # randomly assign zero or one as coefficients to the variables
#' beta <- sample(c(0, 1), 10, replace = TRUE)
#' y <- X %*% beta + rnorm(100)
#'
#' # run the knockoff filter
#' result <- plsko_filter(X, y, q = 0.1)
#' print(result)
#'
#' # compare with the true coefficients
#' which(beta != 0)
#'
#' # Example 2: binary response
#' set.seed(1)
#' X <- matrix(rnorm(100*10), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' # randomly assign zero or one as coefficients to the variables
#' beta <- sample(c(0, 1), 10, replace = TRUE)
#' y <- rbinom(100, 1, plogis(X %*% beta))
#'
#' # run the knockoff filter
#' result <- plsko_filter(X, y, q = 0.1, method = "lasso.logistic")
#' print(result)
#' # compare with the true coefficients
#' which(beta != 0)
#'
#' # Example 3: with neighbourhood information
#' set.seed(1)
#' X <- matrix(rnorm(100*10), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' # randomly assign zero or one as coefficients to the variables
#' beta <- sample(c(0, 1), 10, replace = TRUE)
#' y <- X %*% beta + rnorm(100)
#'
#' # define the neighbourhood information
#' nb.list <- list(1:2, 2:3, 3:4, 4:5, 5:6, 6:7, 7:8, 8:9, 9:10, NULL) # a list of length 10 ((although in this example we know all variables in X are actually independent)
#' # run the knockoff filter
#' result <- plsko_filter(X, y, q = 0.1, nb.list = nb.list)
#' print(result)
#' # compare with the true coefficients
#' which(beta != 0)
#'
#'
plsko_filter <- function(X, y, q = 0.05, method = "lasso.lcd", offset = 0,
                         nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = NULL, sparsity = 1, rmax = 5,
                         seed = 1
                         ){
  set.seed(seed)
  Xk= plsko(X, nb.list = nb.list, threshold.abs = threshold.abs, threshold.q = threshold.q, ncomp = ncomp, sparsity = sparsity, rmax = rmax)
  result = ko_filter(X, Xk, y, q = q, method = method, offset = offset)

  # restructure the result
  result$call <- match.call()
  return(result)
}
