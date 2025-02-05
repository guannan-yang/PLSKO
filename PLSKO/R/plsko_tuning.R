#' @title PLSKO tuning function
#'
#' @import progress
#' @import doParallel
#' @import foreach
#' @import ggplot2
#' @import ggpubr
#' @description Tuning function based on semi-simulation on a series of parameters in PLSKO, including ncomp, threshold.q, threshold.abs and sparsity. The function will return the optimal configuration based on the minimum mean or median FDP over repetitions.
#' @param X The dataset to be used for tuning. A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param p_s The number of assumed important variables to be used to generate artificial y. Default is 10 \% of the number of variables.
#' @param n_ko The number of repetitions for the semi-simulation. Default is 10.
#' @param q The target FDR level. Default is 0.05.
#' @param parallel Logical. If TRUE, run the semi-simulation in parallel. Default is TRUE.
#' @param ncore An integer specifying the number of cores to use for parallel computation. Default is NULL.
#' @param ncomp A integer or a vector of integers to be tested, specifying the number of components to be used in PLSKO. Default is seq(3, 9, 2).
#' @param threshold.q A numeric value or a vector of numeric values (between 0 - 1) to be tested, specifying a quantile threshold to define neighborhoods. Default is 0. This parameter will be used if threshold.abs is not provided.
#' @param threshold.abs A numeric value or a vector of numeric values (between 0 - 1) to be tested, specifying an absolute correlation threshold to define neighborhoods. Default is 0.
#' @param sparsity A numeric value or a vector of numeric values (between 0 (excluded) - 1) to be tested, specifying the sparsity level in the PLS regression. Default is 1.
#' @param simpls Logical. If TRUE, use the SIMPLS algorithm to fit the PLS model. Default is FALSE.
#' @param seed An integer seed for reproducibility. Default is 1.
#' @param fdp.measure A character specifying the measure to be used to select the optimal configuration. Once the mean or median FDP is below the target FDR, the configuration will be return as optimal. Default is "median". Options are "mean", "median", "either" and "both".
#' @param early.stop Logical. If TRUE, the tuning will stop once the mean or median FDP is below the target FDR. Default is FALSE.
#' @return An object of class "plsko.tuning". This object is a list containing at least the following components:
#' \describe{
#'  \item{`full`}{The full result of the tuning.}
#'  \item{`mean`}{The average result of the tuning.}
#'  \item{`median`}{The median result of the tuning.}
#'  \item{`optimal`}{The optimal configuration based on the minimum mean or median FDP over repetitions.
#'  Each component is a dataframe containing the configuration, performance metrics (including FDP and TPP), run times and average performance metrics.}
#'  }
#'
#' @examples
#' data(cfRNA_placenta)
#' X <- cfRNA_placenta$counts
#' cfRNA_tune <- plsko_tuning(X, p_s = 10, n_ko = 10, q = 0.05, parallel = TRUE, ncore = NULL, ncomp = seq(3, 9, 2), threshold.q = 0, threshold.abs = c(0.5, 0.3, 0), sparsity = 1, simpls = F, seed = 1, fdp.measure = "median", early.stop = F)
#' plot(cfRNA_tune)
#'
#' @export
plsko_tuning <- function(X, p_s = round(0.1 * ncol(X)), n_ko = 10, q = 0.05, parallel = TRUE, ncore = NULL,
                         ncomp = seq(3, 9, 2), threshold.q = 0, threshold.abs = 0, sparsity = 1, simpls = F,
                         seed = 1, fdp.measure = "median", early.stop = F) {

  set.seed(seed)
  X <- scale(X)
  if(is.null(colnames(X))) colnames(X) <- as.character(1:ncol(X))
  if(is.null(rownames(X))) rownames(X) <- as.character(1:nrow(X))


  # transfer threshold.q to equivalent threshold.abs
  if(!is.null(threshold.q) & is.null(threshold.abs)){
    cor_mat <- cor(X)
    threshold.abs = quantile(unlist(abs(cor_mat)), prob = threshold.q, names = F)
  }
  else if(is.null(threshold.q) & is.null(threshold.abs)){
    threshold.abs = 0
  }

  # sort the configurations
  ## ncomp test as ascending order
  ncomp <- sort(ncomp)
  ## threshold.abs test as descending order -- the higher the better
  threshold.abs <- sort(threshold.abs, decreasing = T)
  ## no preferable order for sparsity

  test.grid <- expand.grid(ncomp = ncomp, threshold.abs = threshold.abs, sparsity = sparsity)

  # Test across the grid
  full.result <- data.frame()
  mean.result <- data.frame()
  median.result <- data.frame()

  for (i in 1:nrow(test.grid)) {

    ncomp = test.grid$ncomp[i]
    threshold.abs = test.grid$threshold.abs[i]
    sparsity = test.grid$sparsity[i]

    # run semi-simulation
    semi_result <- plsko_semi_sim(X = X, p_s = p_s, q = q, n_ko = n_ko, plsko.ncomp = ncomp, plsko.threshold.abs = threshold.abs, plsko.sparsity = sparsity, simpls = simpls, seed = seed, parallel = parallel, ncores = ncore) # set seed to make sure every configuration is tested on the same y

    # store the result
    full.result <- rbind(full.result, semi_result$res)
    mean.result <- rbind(mean.result, semi_result$average.res)
    median.result <- rbind(median.result, semi_result$median.res)

    if(early.stop){
      if(fdp.measure == "mean"){
        if(semi_result$average.res$mean.fdp < q) {
          print(paste0("Config", i, ":", "ncomp = ", ncomp, ", threshold.abs = ", threshold.abs, ", sparsity = ", sparsity, ", got the mean FDP of ", semi_result$average.res$mean.fdp, "under the target ", q, ". Stop tune."))
          break
        }

      }
      else if(fdp.measure == "median"){
        if(semi_result$median.res$median.fdp < q) {
          print(paste0("Config", i, ":", "ncomp = ", ncomp, ", threshold.abs = ", threshold.abs, ", sparsity = ", sparsity, ", got the median FDP of", semi_result$median.res$median.fdp, "under the target", q, ". Stop tune."))
          break
        }
      }
      else if(fdp.measure == "either"){
        if(semi_result$average.res$mean.fdp < q | semi_result$median.res$median.fdp < q){
          print(paste0("Config", i, ":", "ncomp = ", ncomp, ", threshold.abs = ", threshold.abs, ", sparsity = ", sparsity, ", got the mean FDP of", semi_result$average.res$mean.fdp, "and the median FDP of", median.result$median.fdq, "under the target", q, ". Stop tune"))
          break
        }
        else if(fdp.measure == "both"){
          if(semi_result$average.res$mean.fdp < q & semi_result$median.res$median.fdp < q){
            print(paste0("Config", i, ":", "ncomp = ", ncomp, ", threshold.abs = ", threshold.abs, ", sparsity = ", sparsity, ", got the mean FDP of", semi_result$average.res$mean.fdp, "and the median FDP of", median.result$median.fdq, "under the target", q, ". Stop tune"))
            break
          }
        }
        else{
          stop("fdp.measure should be either 'mean', 'median', 'either' or 'both'")
        }
      }
    }
  }
  # Search for the minimum FDP
  if(fdp.measure == "mean"){
    optimal.result <- mean.result[which.min(mean.result$mean.fdp),]
    print(paste0("Optimal configuration: ncomp = ", optimal.result$ncomp, ", threshold.abs = ", optimal.result$threshold.abs, ", sparsity = ", optimal.result$sparsity, ", mean FDP = ", optimal.result$mean.fdp))
  }
  else if(fdp.measure == "median"){
    optimal.result <- median.result[which.min(median.result$median.fdp),]
    print(paste0("Optimal configuration: ncomp = ", optimal.result$ncomp, ", threshold.abs = ", optimal.result$threshold.abs, ", sparsity = ", optimal.result$sparsity, ", median FDP = ", optimal.result$median.fdp))
  }
  else if(fdp.measure == "either" | fdp.measure == "both"){
    if(mean.result$fdp < median.result$fdp){
      optimal.result <- mean.result[which.min(mean.result$mean.fdp),]
      print(paste0("Optimal configuration: ncomp = ", optimal.result$ncomp, ", threshold.abs = ", optimal.result$threshold.abs, ", sparsity = ", optimal.result$sparsity, ", mean FDP = ", optimal.result$mean.fdp))
    }
    else{
      optimal.result <- median.result[which.min(median.result$median.fdp),]
      print(paste0("Optimal configuration: ncomp = ", optimal.result$ncomp, ", threshold.abs = ", optimal.result$threshold.abs, ", sparsity = ", optimal.result$sparsity, ", median FDP = ", optimal.result$median.fdp))
    }
  }

  full.result$optimal <- F
  full.result$optimal[full.result$ncomp == optimal.result$ncomp & full.result$threshold.abs == optimal.result$threshold.abs & full.result$sparsity == optimal.result$sparsity] <- T

  # object type as plsko.tuning
  result <- list(full = full.result, mean = mean.result, median = median.result, optimal = optimal.result)
  class(result) <- "plsko.tuning"
  return(result)
}

#' @title Plot the tuning result
#' @description Plot the tuning result as a boxplot of FDP and TPP, with the optimal configuration highlighted with a 'star' on the correspondence coordinate in the plot.
#' @param plsko.tuning The tuning result from plsko_tuning function.
#' @param target.fdr The target FDR to be used as a dash line in the plot. Default is 0.05.
#' @export
plot.plsko.tuning <- function(plsko.tuning, target.fdr = 0.05){
  # plot full.result as a boxplot of FDP: x-axis ncomp, group by threshold.abs, panel by sparsity, with target FDR dash line and optimal configuration highlighted with a 'star' on the correspondence coordinate in the plot

  full.result <- plsko.tuning$full

  full.result$ncomp <- as.factor(full.result$ncomp)
  full.result$threshold.abs <- as.factor(full.result$threshold.abs)
  full.result$sparsity <- as.factor(full.result$sparsity)

  fdp.boxplot <- ggplot(full.result, aes(x = ncomp, y = fdp, fill = threshold.abs)) +
    geom_boxplot() +
    facet_wrap(~sparsity, ncol = 1, labeller = "label_both") +
    geom_hline(yintercept = target.fdr, linetype = "dashed", color = "red") +
    stat_summary(aes(group= threshold.abs, colour = threshold.abs), fun="mean", geom="line") +
    geom_point(data = full.result[full.result$optimal,], aes(x = ncomp, y = 0.8), shape = 8, size = 3, color = "purple", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,0.05,0.1,0.2,0.4,0.6,0.8,1), limits = c(0, 1))+
    labs(x = "ncomp", y = "FDP") +
    theme_minimal()

  tpp.boxplot <- ggplot(full.result, aes(x = ncomp, y = tpp, fill = threshold.abs)) +
    geom_boxplot() +
    facet_wrap(~sparsity, ncol = 1, labeller = "label_both") +
    stat_summary(aes(group= threshold.abs, colour = threshold.abs), fun="mean", geom="line") +
    geom_point(data = full.result[full.result$optimal,], aes(x = ncomp, y = 0.8), shape = 8, size = 3, color = "purple", show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,0.1,0.2,0.4,0.6,0.8,1), limits = c(0, 1))+
    labs(x = "ncomp", y = "TPP") +
    theme_minimal()

  ggarrange(fdp.boxplot, tpp.boxplot, ncol = 2, common.legend = T)
}

#' @title Best seed function
#' @description Return the plsko.seed from the optimal configuration based on the either highest TPP when FDP is under the target, or the lowest FDP, or the lowest squared sum of FDP and (1-TPP)
#' @param plsko.tuning The tuning result from plsko_tuning function.
#' @param seed.measure The measure to be used to select the best seed. Default is "TPP".
#' @param target.fdr The target FDR to be used as a dash line in the plot. Default is 0.05.
#' @keywords internal
best.seed <- function(plsko.tuning, seed.measure = "TPP", target.fdr = 0.05){
  # subset of the optimal configuration
  optimal.full <- plsko.tuning$full[plsko.tuning$full$optimal,]

  # if the minimum FDP is above the target
  if(min(optimal.full$fdp) > target.fdr & seed.measure == "TPP"){
    warning("The minimum FDP is above the target FDR. The seed will be selected based on sum of the squared FDP and (1 - TPP)")
    seed.measure = "FDP_TPP"
  }

  if(seed.measure == "TPP"){
    optimal.seed <- optimal.full[optimal.full$fdp < target.fdr,]$plsko.seed[which(optimal.full$tpp == max(optimal.full$tpp))]
  }
  else if(seed.measure == "FDP"){
    optimal.seed <- optimal.full$plsko.seed[which(optimal.full$fdp == min(optimal.full$fdp))]
  }
  else if(seed.measure == "FDP_TPP"){
    square_sum <- optimal.full$fdp^2 + (1-optimal.full$tpp)^2
    optimal.seed <- optimal.full$plsko.seed[which(square_sum == min(square_sum))]
  }
  else{
    stop("seed.measure should be either 'TPP', 'FDP' or 'FDP_TPP'")
  }
  return(optimal.seed)
}

#' @title Semi-simulation for PLSKO
#' @description Return FDP and TPP of semi-simulation for n_ko repetitions at this configuration
#' @param X The dataset to be used for tuning. A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param p_s The number of assumed important variables to be used to generate artifical . Default is 10 \% of the number of variables.
#' @param q The target FDR to be used as a dash line in the plot. Default is 0.05.
#' @param n_ko The number of repetitions for the semi-simulation. Default is 10.
#' @param plsko.ncomp The number of components to be used in PLSKO.
#' @param plsko.threshold.abs The threshold for the absolute value of the correlation between the original and knockoff variables.
#' @param plsko.sparsity The sparsity level in the PLS regression.
#' @param simpls Logical. If TRUE, use the SIMPLS algorithm to fit the PLS model. Default is FALSE.
#' @param seed An integer seed for reproducibility. Default is 1.
#' @param parallel Logical. If TRUE, run the semi-simulation in parallel. Default is TRUE.
#' @param ncores An integer specifying the number of cores to use for parallel computation. Default is NULL.
#'
#' @return A list containing the result as a dataframe including the configuration, performance metrics (including FDP and TPP), run times and average performance metrics
#'
#' @export
plsko_semi_sim <- function(X, p_s, q, n_ko, plsko.ncomp, plsko.threshold.abs, plsko.sparsity, simpls,
                           seed = Sys.time(), parallel = TRUE, ncores = NULL){

  if(parallel){
    if (!requireNamespace('doParallel', quietly=T)) {
      warning('doParallel is not installed. Without parallelisation, the multiple knockoff sets will be slower to generate', call.=F,immediate.=T)
      parallel=F
    }
    if (!requireNamespace('foreach', quietly=T)) {
      warning('foreach is not installed. Without parallelisation, the multiple knockoff sets will be slower to generate', call.=F,immediate.=T)
      parallel=F
    }

    # Register cores for parallel computation
    if (parallel) {
      all_cores = parallel::detectCores(all.tests = TRUE, logical = TRUE)-2
      if(is.null(ncores)) ncores = all_cores # if not specified, use all cores except one
      if (ncores > all_cores ) {
        warning(paste("The requested number of cores is not available. Using instead",all_cores,"cores"),immediate.=T)
        ncores = all_cores
      }
      if (ncores>1) {
        doParallel::registerDoParallel(cores=ncores)
        parallel = TRUE
      }
      else {
        parallel = FALSE
      }
    }
  }

  # run knockoff filter parallel or not
  if(parallel){
    para.result <- foreach::foreach(i = 1:n_ko, .packages = c('knockoff', 'progress')) %dopar% {

      set.seed(seed+i-1)

      Gy <- DGPy.AR1(X = X, s = p_s, c = 0, y.dis = "Normal", A = c(3, 5))
      y <- Gy$y

      # Generate PLSKO knockoff
      start.time <- Sys.time()
      ko = plsko(X, seed = seed+i-1, ncomp = plsko.ncomp, threshold.abs = plsko.threshold.abs, sparsity = plsko.sparsity, simpls = simpls)
      end.time <- Sys.time()
      run.time <- as.numeric(difftime(end.time, start.time, units = "secs"))

      # Calculate the apply knockoff filter
      S <- ko_filter(X = X, Xk = ko, y = y, q = q, w.method = "lasso.lcd", offset = "both", cores = 1)

      # Calculate the FDP and TPP for this ko
      fdp <- fdr(S$selected, Gy$Beta)
      fdp.plus <- fdr(S$selected.plus, Gy$Beta)
      tpp <- pow(S$selected, Gy$Beta)
      tpp.plus <- pow(S$selected.plus, Gy$Beta)
      plsko.seed <- seed+i-1

      # release memory
      rm(ko)

      res <- list(fdp = fdp, tpp = tpp, fdp.plus = fdp.plus, tpp.plus = tpp.plus, run.time = run.time, plsko.seed = plsko.seed)

    }
    doParallel::stopImplicitCluster()

    fdp = unlist(lapply(para.result, function(x) x$fdp))
    tpp = unlist(lapply(para.result, function(x) x$tpp))
    fdp.plus = unlist(lapply(para.result, function(x) x$fdp.plus))
    tpp.plus = unlist(lapply(para.result, function(x) x$tpp.plus))
    run.time = unlist(lapply(para.result, function(x) x$run.time))
    plsko.seed = unlist(lapply(para.result, function(x) x$plsko.seed))
  }
  else {
    fdp = numeric(n_ko)
    fdp.plus = numeric(n_ko)
    tpp = numeric(n_ko)
    tpp.plus = numeric(n_ko)
    run.time = numeric(n_ko)
    plsko.seed = numeric(n_ko)

    for (i in 1:n_ko) {
      set.seed(seed + i-1 )
      Gy <- DGPy.AR1(X = X, s = p_s, c = 0, y.dis = "Normal", A = c(3, 5))
      y <- Gy$y

      # Generate PLSKO knockoff
      start.time <- Sys.time()

      ko = plsko(X, seed = seed+i-1, ncomp = plsko.ncomp, threshold.abs = plsko.threshold.abs, sparsity = plsko.sparsity, simpls = simpls)

      end.time <- Sys.time()

      # Calculate the apply knokcoff filter
      S <- ko_filter(X = X, Xk = ko, y = y, q = q, w.method = "lasso.lcd", offset = "both", cores = 1)

      # free memory
      rm(ko)

      fdp[i] <- fdr(S$selected, Gy$Beta)
      fdp.plus[i] <- fdr(S$selected.plus, Gy$Beta)
      tpp[i] <- pow(S$selected, Gy$Beta)
      tpp.plus[i] <- pow(S$selected.plus, Gy$Beta)
      run.time[i] <- as.numeric(difftime(end.time, start.time, units = "secs"))
      plsko.seed[i] <- seed + i-1
    }
  }

  # return the result as a dataframe including the configuration, performance metrics, run times and average performance metrics
  res <- data.frame(plsko.seed = plsko.seed,
                    ncomp = plsko.ncomp,  threshold.abs = plsko.threshold.abs, sparsity = plsko.sparsity,
                    fdp = fdp, tpp = tpp, fdp.plus = fdp.plus, tpp.plus = tpp.plus, run.time = run.time)

  average.res <- data.frame(ncomp = plsko.ncomp,threshold.abs = plsko.threshold.abs, sparsity = plsko.sparsity, mean.fdp = mean(fdp), mean.tpp = mean(tpp), mean.fdp.plus = mean(fdp.plus), mean.tpp.plus = mean(tpp.plus), mean.run.time = mean(run.time))
  median.res <- data.frame(ncomp = plsko.ncomp, threshold.abs = plsko.threshold.abs, sparsity = plsko.sparsity, median.fdp = median(fdp), median.tpp = median(tpp), median.fdp.plus = median(fdp.plus), median.tpp.plus = median(tpp.plus), median.run.time = median(run.time))

  result <- list(res = res, average.res = average.res, median.res = median.res)
  return(result)
}

#' @title Calculate False discovery proportion
#' @description Calculate the false discovery proportion (FDP) based on the selected variables and the true support.
#' @param S The selected variables.
#' @param beta.true The true support.
#' @return The false discovery proportion (FDP).
#' @keywords internal
fdr <- function(S, beta.true) {
  fdp = sum(beta.true[S] == 0)/max(1, length(S))
  return(fdp)
}

#' @title Calculate True positive proportion
#' @description Calculate the true positive proportion (TPP) based on the selected variables and the true support.
#' @param S The selected variables.
#' @param beta.true The true support.
#' @return The true positive proportion (TPP).
#' @keywords internal
pow <- function(S, beta.true) {
  tpp = sum(beta.true[S] != 0)/sum(beta.true != 0)
  return(tpp)
}

#' @title Generate artificial response variable y
#'
#' @description Generate artificial response variable y based on the true support and coefficients.
#'
#' @param X The dataset to be used for tuning. A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param s The number of assumed important variables to be used to generate artifical .
#' @param A The range of the coefficients for the true support.
#' @param c The noise level. Default is 0.
#' @param y.dis The distribution of the response variable. Default is "Normal".
#'
#' @return A list containing the true support and the generated response variable y.
#'
#' @keywords internal
DGPy.AR1 <- function(X,s,A, c = 0, y.dis = "Normal"){
  n = nrow(X)
  p = ncol(X)
  B.true <- rep(0,p)
  S.true <- sample(1:p, size = s)
  S.true <- sort(S.true) # True Support
  C.unif <- c(runif(p, A[1], A[2]), runif(p, -A[2], -A[1]))
  C <- sample(C.unif, size = s , replace = TRUE) # coefficients for true signals
  B.true [S.true] <- C # A true vector of coefficients
  epsilon <- matrix(rnorm(n),n,1)  # Error Term

  linear.fit <- X %*% B.true # Response variable

  if(y.dis == "Normal"){
    y <- linear.fit + sqrt(c*s)*epsilon
    # noise.ratio = var(sqrt(c*s)*epsilon)/var(linear.fit)
    # y = list(y, noise.ratio)
  }
  else if(y.dis == "NB"){
    mu = 2^linear.fit*400
    prob = 0.97
    size = prob/(1-prob)*mu
    y <- rnbinom(n, size = size, mu = mu)
  }

  else if(y.dis == "Binary"){
    p.logit <- 1/(1+exp(-(linear.fit-mean(linear.fit))))# logistics
    y <- rbinom(n, size = 1, prob = p.logit)
  }

  obj = list(Beta = B.true, y = y)
  # }
  return(obj)
}
