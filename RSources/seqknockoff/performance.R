#' False discovery proportion (fdp) as function of selection and known negatives:
#'
#' @param selected vector of indices of selected variables
#' @param negatives vector of indices of known null variables (that don't influence response)
#'
#' @return false discovery rate
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' eval_fdp(selected=c(1,2,3,5), negatives=5:10)
eval_fdp <- function(selected, negatives) {

  fdr <- ifelse(length(selected) > 0, sum(selected %in% negatives)/length(selected), 0)

  return(fdr)

}

#' True positive proportion (tpp) as function of selection and known positives:
#'
#' @param selected vector of indices of selected variables
#' @param negatives vector of indices of known non-null variables (that influence response)
#'
#' @return true positive rate
#' @export
#'
#' @examples
#' #' library(seqknockoff)
#'
#' eval_tpp(selected=c(1,2,3,5), positives=1:4)
eval_tpp <- function(selected, positives) {

  tpr <- sum(positives %in% selected)/length(positives)

  return(tpr)

}
