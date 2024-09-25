
#' Knockoff filter with generated knockoff variable set
#'
#' @import knockoff
#' @import glmnet
#'
#' @param X A \eqn{n \times p} numeric matrix of predictors.
#' @param Xk A \eqn{n \times p} numeric matrix of knockoff predictors.
#' @param y A numeric or factor type repose vector of length \eqn{n}.
#' @param q A numeric value of the target false discovery rate (FDR) level.
#' @param w.method A character string specifying the method to compute feature importance statistics. Default is \code{"lasso.lcd"}. Other options include \code{"lasso.logistic"} for binary response variable, \code{"lasso.max.lambda"} for the maximum lambda value for the first entry on the path, and \code{"RF"} for random forest. See \code{\link{ko_filter}} or \code{\link{knockoff::knockoff.filter}} for more details.
#' @param offset A numeric or string value to adjust the knockoff threshold. Default is \eqn{0} (control modified FDR). Other options include \eqn{1} (yields a slightly more conservative procedure ("knockoffs+") that controls the FDR according to the usual definition) and \code{"both"} (returns results from both "knockoffs" and "knockoffs+").
#' @param ... Additional arguments passed to the underlying methods. See functions `stat.lasso_coefdiff`, `stat.lasso_lambdadiff`, `stat.lasso_coefdiff_bin`, `stat.random_forest` in `knockoff` package for more details. More options will be provided in the future.
#'
#' @return An object of class "knockoff.result". This object is a list
#'  containing at least the following components:
#' \describe{
#'   \item{`call`}{The matched call.}
#'   \item{`statistic`}{The computed test statistics `W` of length \eqn{p}.}
#'   \item{`threshold`}{The knockoff threshold used to determine discoveries.}
#'   \item{`selected`}{The named vector of selected variables.}
#' }
#'
#' @references
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://web.stanford.edu/group/candes/knockoffs/index.html}{https://web.stanford.edu/group/candes/knockoffs/index.html}
#'
#'   Barber and Candes,
#'   Controlling the false discovery rate via knockoffs.
#'   Ann. Statist. 43 (2015), no. 5, 2055--2085.
#'
#' @examples
#' # Example usage of ko.filter
#' set.seed(123)
#' X <- matrix(rnorm(100*10), 100, 10)
#' colnames(X) <- paste0("X", 1:10)
#' Xk <- plsko(X) # generate knockoff variables by PLSKO
#'
#' # Example 1: continuous response
#' # randomly assign zero or one as coefficients to the variables
#' beta <- sample(c(0, 1), 10, replace = TRUE)
#' y <- X %*% beta + rnorm(100)
#'
#' # run the knockoff filter
#' result <- ko_filter(X, Xk, y, q = 0.1)
#' print(result)
#'
#' # compare with the true coefficients
#' which(beta != 0)
#'
#' # Example 2: knockoff+
#' result.plus <- ko_filter(X, Xk, y, q = 0.1, offset = 1)
#' print(result.plus) # return NULL since no variables are selected as a conservative procedure
#'
#' # Example 3: binary response
#' y.bin <- rbinom(100, 1, 1/(1+ exp(-(X %*% beta)))) # convert to binary response
#' result.bin <- ko_filter(X, Xk, y.bin, q = 0.1,w.method = "lasso.logistic")
#' print(result.bin)
#'
#' @export
ko_filter <- function(X, Xk, y, q = 0.05,w.method = "lasso.lcd", offset = 0, ...){

  n = nrow(X)
  p = ncol(X)
  if (is.data.frame(X)) {
    X.names = colnames(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')

  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if( is.numeric(y) ) y = as.vector(y)

  # dimension check
  if(nrow(X) != length(y)){
    stop('The number of rows in X must be equal to the length of y')
  }
  if(ncol(X) != ncol(Xk)){
    stop('The number of columns in X must be equal to the number of columns in Xk')
  }
  if(nrow(X) != nrow(Xk)){
    stop('The number of rows in X must be equal to the number of rows in Xk')
  }

  if(offset!=1 && offset!=0 && offset!="both") {
    stop('Input offset must be either 0, 1, or both')
  }

  if(w.method == "lasso.lcd" & is.factor(y)){
    warning('lasso.lcd is not applicable for categorical response variable, change to lasso.logistic')
    w.method = "lasso.logistic"
  }
  if(w.method == "lasso.logistic" & !is.factor(y)){
    warning('lasso.logistic is not applicable for continuous response variable, change to lasso.lcd')
    w.method = "lasso.lcd"
  }

  if(w.method == "lasso.lcd"){
    W <- stat.lasso_coefdiff(X, Xk, y, ...)
  }
  else if(w.method == "lasso.max.lambda"){
    W <- stat.lasso_lambdadiff(X, Xk, y, ...)
  }
  # else if(w.method == "pls.regression.lcd"){
  #   W <- W_spls_lcd(X, Xk, y, lcd_comp = 1)
  # }
  #
  # else if(w.method == "plsda.lcd"){ # exclude when continuous y design
  #   W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1)
  # }
  # else if(w.method == "spls.regression.lcd"){
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_spls_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }
  # else if(w.method == "splsda.lcd"){ # exclude when continuous y design
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_splsda_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }

  else if(w.method == "lasso.logistic"){ # exclude when continuous y design
    W <- stat.lasso_coefdiff_bin(X, Xk, y, ...)
  }

  else if(w.method == "RF"){
    W <- stat.random_forest(X, Xk, y)
  }

  result <- ko_withW(W, q = q, offset = offset, X.names = X.names)

  # return the function call
  result$call <- match.call()
  return(result)
}

#' Knockoff Filtering with Precomputed Test Statistics
#'
#' Computes the knockoff threshold and selects variables based on the computed test statistics.
#'
#' @param W A (named) numeric vector of test statistics. The length of the vector should be equal to the number of variables. The output will be names if this vector is named.
#' @param q A numeric value specifying the false discovery rate (FDR) level. Default is `0.05`.
#' @param offset  A numeric or string value to adjust the knockoff threshold. Default is \eqn{0} (control modified FDR). Other options include \eqn{1} (yields a slightly more conservative procedure ("knockoffs+") that controls the FDR according to the usual definition) and \code{"both"} (returns results from both "knockoffs" and "knockoffs+").
#' @param X.names Optional. A character vector of variable names. If provided, the output will be named.
#'
#' @return An object of class "knockoff.result". This object is a list
#'  containing at least the following components:
#' \describe{
#'   \item{`call`}{The matched call.}
#'   \item{`statistic`}{The computed test statistics `W` of length \eqn{p}.}
#'   \item{`threshold`}{The knockoff threshold used to determine discoveries.}
#'   \item{`selected`}{The named vector of selected variables.}
#' }
#'
#' @examples
#' # Example showing self-defined test statistics
#' set.seed(123)
#' n = 100
#' p = 10
#' # generate random test data that in one group, some random selected variables have a higher mean
#' X = matrix(rnorm(n*p), n, p)
#' colnames(X) = paste0("X", 1:p)
#' y = rep(0:1, each = n/2) #assign group
#' beta = sample(c(0, 1), p, replace = TRUE) # randomly assign zero or one as which variables have higher mean
#' X[y==1, beta==1] = X[y==1, beta==1] + rnorm(n/2, mean = 2) # add 2 to the mean of the selected variables
#'
#' Xk <- plsko(X)
#'
#' # generate self-defined test statistics
#' Z <- abs(apply(X[y==1,], 2, mean) - apply(X[y==0,], 2, mean)) # use difference of means as test statistics
#' Zk <- abs(apply(Xk[y==1,], 2, mean) - apply(Xk[y==0,], 2, mean))
#' W <- Z - Zk # compute the difference of test statistics between original and knockoff variables
#' names(W) <- colnames(X)
#'
#' # run the knockoff filter
#' result <- ko_withW(W, q = 0.1)
#' print(result)
#' which(beta != 0) # compare with the true coefficients
#'
#'
#' @export
ko_withW <- function(W, q = 0.05, offset = 0, X.names = NULL){
  if(!is.numeric(W)) stop('Input W must be a numeric vector')
  if(!is.null(X.names) && length(X.names) != length(W)){
    stop('Length of X.names must be the same as the length of W')
  }

  if(!is.null(names(W))& is.null(X.names)) X.names = names(W)

  p = length(W)
  # find the knockoff threshold T
  t = sort(c(0, abs(W)))

  ratio = c(rep(0, p))
  ratio.plus = c(rep(0, p))

  for (j in 1:p) {
    if(offset==0){
      ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }
    else if(offset == 1){
      ratio[j] = (1+sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }
    else if(offset == "both"){
      ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
      ratio.plus[j] = (1+sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }
  }

  id = which(ratio <= q)[1]

  if(length(id) == 0){
    T = Inf
  } else {
    T = t[id]
  }
  # set of discovered variables
  S = integer0_test(which(W >= T))

  if (!is.null(X.names) & !is.null(S))
    names(S) = X.names[S]

  if(offset == "both"){
    id.plus = which(ratio.plus <= q)[1]
    if(length(id.plus) == 0){
      T.plus = Inf
    } else {
      T.plus = t[id.plus]
    }
    S.plus = integer0_test(which(W >= T.plus))
    if (!is.null(X.names))
      names(S.plus) = X.names[S.plus]
  }

  result <- structure(list(call = match.call(),
                           #X = X,
                           #Xk = Xk,
                           #y = y,
                           statistic = W,
                           threshold = T,
                           selected = S),
                      class = 'knockoff.result')
  if(offset == "both"){
    structure(list(call = match.call(),
                   statistic = W,
                   threshold = T,
                   threshold.plus = T.plus,
                   selected = S,
                   selected.plus = S.plus),
              class = 'knockoff.result')
  }

  return(result)
}

#' Print results for the knockoff filter
#'
#' Prints the list of variables selected by the knockoff filter and the corresponding function call.
#'
#' @param x the output of a call to \code{ko_filter}
#' @param ... unused
#'
#' @method print knockoff.result
#' @export
print.knockoff.result <- function(x, plus = F, ...) {
  cat('Call:\n')
  print(x$call)
  cat('\nSelected variables:\n')
  print(x$selected)

  if(plus){
    cat('\nSelected variables (knockoffs+):\n')
    print(x$selected.plus)
  }
}


#' Integer Zero Test
#'
#' Helper function to handle empty selections.
#'
#' @param which A vector of indices.
#'
#' @return Either the input indices or `NULL` if the input is empty.
#' @keywords internal

integer0_test <- function(which){
  if(identical(which, integer(0))){
    return(NULL)
  }
  else {return(which)}
}
