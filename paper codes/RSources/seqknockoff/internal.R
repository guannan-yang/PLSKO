#' Internal check of whether class of response is compatible with family
#'
#' Do not call this function on its own
#'
#' @param y response variable
#' @param family character string, either "gaussian" or "binomial"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @examples
check_family <- function(y, family) {

  if (!((is.numeric(y) & family=="gaussian") | (is.factor(y) & length(unique(y))==2 & family=="binomial"))) {
    stop("Either 1) y should be numeric and family='gaussian', or 2) y should be a binary factor and family='binomial'")
  }

}

#' Internal check of whether inpu data frame (or tibble) is of the right format
#'
#' Do not call this function on its own
#'
#' @param X data frame or tibble
#' @param method character string, either "seq" or "mx"
#'
#' @return this function doesn't return anything
#' @export
#'
#' @examples
check_design <- function(X, method="seq") {

  if(!("data.frame" %in% class(X))) {
    stop("X should be either a data.frame or tibble")
  }

  if(method=="seq" & ncol(X)<=2) {
    stop("X should have ncol(X) > 2")
  }

  if(method=="seq" & sum(!unlist(lapply(X, class)) %in% c("factor", "numeric")) > 0) {
    stop("X should only contain columns of class 'numeric' or 'factor'")
  }

  if(method=="mx" & sum(!unlist(lapply(X, class)) %in% c("numeric")) > 0) {
    stop("X should only contain columns of class 'numeric'")
  }

}



#' Select variables based on (heuristic) mode of multiple variable selections
#'
#' Do not call this function on its own
#'
#' @param S list of variable selection indices
#' @param p number of variables. Each element of the list of selection indices should be a subset of 1:p.
#' @param trim trimming probability threshold. A sensible default is \code{trim=0.5}.
#'
#' @return a single "most frequent" variable selection among the multiple selections in S.
#' @export
#'
#' @examples
find_single_optimal_variable_set <- function(S, p, trim=0.5) {

  candidates <- 1:p

  countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

  Nknockoff <- length(S)

  var.freq.table <- table(factor(unlist(S), levels=candidates))

  which.remove <- as.numeric(which(var.freq.table < trim*Nknockoff))

  trimmed.selected <- lapply(S, function(x) paste(setdiff(x, which.remove), collapse=" "))

  model.freq.table <- table(unlist(trimmed.selected))

  best.trimmed.selected <- names(which(model.freq.table==max(model.freq.table)))

  # Resolve ties by choosing most parsimonious model:
  best.single.selected <- best.trimmed.selected[which.min(countSpaces(best.trimmed.selected))]

  selected <- as.integer(unlist(strsplit(best.single.selected, " ")))

  return(selected)

}
