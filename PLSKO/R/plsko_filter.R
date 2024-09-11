
plsko_filter <- function(X, y, q = 0.05, method = "lasso.lcd", offset = 0,
                         nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = NULL, sparsity = 1,
                         seed = 1
                         ){
  set.seed(seed)
  Xk= plsko(X, nb.list = nb.list, threshold.abs = threshold.abs, threshold.q = threshold.q, ncomp = ncomp, sparsity = sparsity)
  result = ko_filter(X, Xk, y, q = q, method = method, offset = offset)

  return(result)
}
