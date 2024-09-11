
#' Knockoff filter with generated knockoff variable set
#'
#' @import knockoff
#' @import glmnet
#'
#' @param X A \eqn{n \times p} numeric matrix of predictors.
#' @param Xk A \eqn{n \times p} numeric matrix of knockoff predictors.
#' @param y A numeric or factor type repose vector of length \eqn{n}.
#' @param q A numeric value of the target false discovery rate (FDR) level.
#' @param method A string specifying the method to compute test statistics. Options include
#'   `"lasso.lcd"`, `"lasso.max.lambda"`, `"lasso.logistic"`, `"RF"`. Default is `"lasso.lcd"`.
#' @param offset A numeric or string value to adjust the knockoff threshold. Default is `0` (control modifed FDR). Other options include `1` (yields a slightly more conservative procedure ("knockoffs+") that controls the FDR according to the usual definition) and `"both"` (returns results from both "knockoffs" and "knockoffs+").
#' @param ... Additional arguments passed to the underlying methods. See functions `stat.lasso_coefdiff`, `stat.lasso_lambdadiff`, `stat.lasso_coefdiff_bin`, `stat.random_forest` in `knockoff` package for more details.
#'
#' @return An object of class "knockoff.result". This object is a list
#'  containing at least the following components:
#' \describe{
#'   \item{`call`}{The matched call.}
#'   \item{`X`}{The original matrix of predictors.}
#'   \item{`Xk`}{The knockoff matrix of predictors.}
#'   \item{`y`}{The response vector.}
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
#' result <- ko.filter(X, Xk, y, q = 0.1)
#'
#' @export
ko.filter <- function(X, Xk, y, q = 0.05, method = "lasso.lcd", offset = 0, ...){

  n = nrow(X)
  p = ncol(X)
  if (is.data.frame(X)) {
    X.names = names(X)
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

  if(offset!=1 && offset!=0 && offset!="both") {
    stop('Input offset must be either 0, 1, or both')
  }


  if(method == "lasso.lcd"){
    W <- stat.lasso_coefdiff(X, Xk, y, ...)
  }
  else if(method == "lasso.max.lambda"){
    W <- stat.lasso_lambdadiff(X, Xk, y, ...)
  }
  # else if(method == "pls.regression.lcd"){
  #   W <- W_spls_lcd(X, Xk, y, lcd_comp = 1)
  # }
  #
  # else if(method == "plsda.lcd"){ # exclude when continuous y design
  #   W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1)
  # }
  # else if(method == "spls.regression.lcd"){
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_spls_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }
  # else if(method == "splsda.lcd"){ # exclude when continuous y design
  #   keepX <- round(ncol(X)*sparsity)+1
  #   W <- W_splsda_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  # }

  else if(method == "lasso.logistic"){ # exclude when continuous y design
    W <- stat.lasso_coefdiff_bin(X, Xk, y, ...)
  }

  else if(method == "RF"){
    W <- stat.random_forest(X, Xk, y)
  }

  result <- ko.withW(W, q = q, offset = offset)

  return(result)
}

#' Knockoff Filtering with Precomputed Test Statistics
#'
#' Computes the knockoff threshold and selects variables based on the computed test statistics.
#'
#' @param W A numeric vector of test statistics.
#' @param q A numeric value specifying the false discovery rate (FDR) level. Default is `0.05`.
#' @param offset A numeric or string value to adjust the knockoff threshold. Default is `0` (control modifed FDR). Other options include `1` (yields a slightly more conservative procedure ("knockoffs+") that controls the FDR according to the usual definition) and `"both"` (returns results from both "knockoffs" and "knockoffs+").
#'
#' @return An object of class "knockoff.result". This object is a list
#'  containing at least the following components:
#' \describe{
#'   \item{`call`}{The matched call.}
#'   \item{`X`}{The original matrix of predictors.}
#'   \item{`Xk`}{The knockoff matrix of predictors.}
#'   \item{`y`}{The response vector.}
#'   \item{`statistic`}{The computed test statistics `W` of length \eqn{p}.}
#'   \item{`threshold`}{The knockoff threshold used to determine discoveries.}
#'   \item{`selected`}{The named vector of selected variables.}
#' }
#'
#' @export
ko.withW <- function(W, q = 0.05, offset = 0){
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

  if (!is.null(X.names))
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
#' @param x the output of a call to ko.filter
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
