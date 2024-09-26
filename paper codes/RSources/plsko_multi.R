# multiple plsko
library(mixOmics)

plsko.multi <- function(X, threshold.abs = NULL, threshold.q = NULL, ncomp = 2, sparsity = 1, gauss.err = F, M.ko = 1){
  start.time <- Sys.time()
  n <- nrow(X)
  p <- ncol(X)
  
  #random swap column of X
  sample.order.mat <- diag(1, p, p)[,sample(p)]
  X <- X %*% sample.order.mat
  
  mu <- colMeans(X)  
  call <- list(thres.abs = threshold.abs,
               thres.q = threshold.q,
               ncomp = ncomp,
               sparsity = sparsity,
               gauss.err = gauss.err)
  
  X <- scale(X, center = T, scale = F)
  
  # Calculate the correlation matrix
  cor_mat <- cor(X)
  # Replace the diagonal values with 0, since they are always 1
  diag(cor_mat) <- 0
  
  if(!is.null(threshold.q)){
    threshold <- quantile(unlist(abs(cor_mat)), prob = threshold.q, names = F)
  } 
  else if(!is.null(threshold.abs)){
    threshold <- threshold.abs
  }
  else {
    threshold <- quantile(unlist(abs(cor_mat)), prob = 0.9, names = F)
  }
  
  neighborhoods <- vector(mode = "list", length = p)
  for (i in 1:p) {
    # Find the indices of the columns with correlation greater than the threshold
    cols <- which(abs(cor_mat[i,]) > threshold)
    
    # Add the indices to the list of neighborhoods
    if (length(cols) > 0) {
      neighborhoods[[i]] <- cols
    }
  }
  
  # if neighbours more than or equal to n, using pls-regression to estimate variable i's distribution;
  # if neighbours less than n, using linear regression
  # if no neighbours, treat as random gaussian then sampling

  X_k.list = list()
  
  for (m in 1:M.ko){
    print(m)
    X <- as.data.frame(X)
    X_k <- matrix(NA,nrow(X),ncol(X))
    X_k <- as.data.frame(X_k)
    rownames(X_k) <- rownames(X)
    X_k <- X
    colnames(X_k) <- paste0(colnames(X_k),"k","_",m)
    
    for (i in 1:p){
      tryCatch({  
        #print(i)
        nb <- neighborhoods[[i]]
        k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
        X.nb <- X[,nb]
        X_k.nb <- X_k[, k.nb]
        
        if(m > 1){# previous knockoff's k.nb
        
        X_k_1_m <- lapply(X_k.list, function(xk){
          xk_nb <- xk[, k.nb]
          return(xk_nb)
        })
        X_k_1_m <- do.call("cbind", X_k_1_m)
        X_k_1_m <- as.data.frame(X_k_1_m)
        if(length(X_k_1_m)!=0){
          names(X_k_1_m) <- paste0(m,"k_",1:ncol(X_k_1_m))
        }
        
        X_k.nb <- cbind(X_k.nb, X_k_1_m)}
        
        if(length(nb) ==0){
          #print("0")
          Y <- X[,i]
          Y.hat <- 0
        }
        
        else if(length(nb)== 1 & (length(nb)*m < n)){
          #print("linear")
          X.run <- cbind(X.nb, X_k.nb)
          Y <- X[,i]
          
          Y.hat <- linear.regression.generator(Y, X.run)
        }
        
        else{
          #if(length(k.nb) < n/2){
          #print("pls")
          X.run <- cbind(X.nb, X_k.nb)
          Y <- X[,i]
          if(is.null(ncomp)){
            this.ncomp <- ceiling(min(n/2, ncol(X.nb)/2))
          }
          else this.ncomp <- ncomp
          
          this.ncomp <- min(this.ncomp, ncol(X.run)-1)
          
          keepX <- rep(round(sparsity*ncol(X.run)), this.ncomp)
          Y.hat <- plsda.recovery.generator(Y, X.run, ncomp = this.ncomp, keepX = keepX)
        }
        
        Y.res <- Y - Y.hat
        # resSD <- sqrt(sum(Y.res^2))
        # res <- rnorm(n, sd = resSD)
        # X_ki <- Y.hat + resSD
        if(!gauss.err){
          # permuation
          res <- sample(Y.res) 
        }
        else {
          sd.err <- sd(Y.res)
          res <- rnorm(n, 0, sd.err)
        }
        
        X_k[,i] <- Y.hat + res
      }
      , error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
    }
    
    X_k <- apply(X_k, 1, function(x){x+mu})
    X_k <- t(X_k)
    
    #swap back
    X_k <- X_k %*% t(sample.order.mat)
    colnames(X_k) <- paste0(colnames(X),"k","_",m)
    
    X_k.list[[m]] <- X_k
  }

  end.time <- Sys.time()
  run.time <- difftime(end.time, start.time, units = "secs")
  obj <- list(X_k.list = X_k.list, 
              call = call,
              run.time = run.time)
  return(obj)  
}

plsda.recovery.generator <- function(Y, X, ncomp, keepX = rep(ncol(X))){
  X <- as.data.frame(X)
  #X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")
  Y.hat <- predict(pls, X)$predict[1:n, 1, ncomp]
  return(Y.hat)
}

linear.regression.generator <- function(Y, X){
  X <- as.data.frame(X)
  X <- X[!duplicated(as.list(X))]
  X <- as.matrix(X)
  betas <- ginv(crossprod(X))%*% t(X) %*% Y
  Y.hat <- X %*% betas
  return(Y.hat)
}

