#source('/data/gpfs/projects/punim0613/Guan/Sim2/RSource/KOBT/create.pc.knockoff.R')
#source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/plsko_plsonly.R")
# source('./RSource/KOBT/create.pc.knockoff.R')
# source("./RSource/plsko_plsonly.R")

hd.ko <- function(X, threshold.abs = NULL, threshold.q = NULL, fit.method = "ginv", r.fit = 3, gauss.err = F){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #random swap column of X
  sample.order.mat <- diag(1, p, p)[,sample(p)]
  X <- X %*% sample.order.mat
  
  mu <- colMeans(X)
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
    threshold <- quantile(unlist(abs(cor_mat)), prob = 0.95, names = F)
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
  
  # if neighbours more than or equal to n, using plsda to estimate variable i's distribution;
  # if neighbours less than n, using linear regression
  # if no neighbours, treat as random gaussian then sampling
  X <- as.data.frame(X)
  X_k <- matrix(NA,nrow(X),ncol(X))
  X_k <- as.data.frame(X_k)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(colnames(X),"k")
  
  for (i in 1:p){
    #tryCatch({  
      #print(i)
      nb <- neighborhoods[[i]]
      k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
      X.nb <- X[,nb]
      X_k.nb <- X_k[, k.nb]
      
      X.run <- cbind(X.nb, X_k.nb)
      Y <- X[,i]
      
      if(length(nb) ==0){
        print("0")
        Y <- X[,1]
        Y.hat <- 0
      }
      
      else{
      if(fit.method == "ginv"|length(nb)==1){
        #X.run <- as.matrix(X)
        #print(eigen(crossprod(as.matrix(X.run)),only.values = T)$values)
        Y.hat <- ginv.regression.generator(Y, X.run)
      }
      else if(fit.method == "pca"){
        Y.hat <- pca.regression.generator(Y, X.run, ncomp = r.fit)
      }
      else if(fit.method == "pls"){
        Y.hat <- plsda.recovery.generator(Y, X.run, ncomp = r.fit)
      }
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
    #, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  #}
  
  X_k <- apply(X_k, 1, function(x){x+mu})
  X_k <- t(X_k)
  
  #swap back
  X_k <- X_k %*% t(sample.order.mat)
  
  obj <- X_k
  return(obj)  
}

ginv.regression.generator <- function(Y, X){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  #betas <- MASS::ginv(crossprod(X))%*% t(X) %*% Y
  betas <- ginv_for_crossprod(X) %*% t(X) %*% Y
  Y.hat <- X %*% betas
  Y.hat[is.na(Y.hat)] <- 0
  Y.hat[is.infinite((Y.hat))] <- 0
  Y.hat[is.null(Y.hat)] <- 0
  return(Y.hat)
}

ginv_for_crossprod <- function(X, tol = sqrt(.Machine$double.eps)){ # gain the ginv(crossprod(X))
  Xsvd <- svd(X)
  if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) Xsvd$v %*% (1/(Xsvd$d)^2 * t(Xsvd$v))
  else if(!any(Positive)) array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop=FALSE] %*% ((1/(Xsvd$d[Positive])^2) * t(Xsvd$v[, Positive, drop=FALSE]))
}

pca.regression.generator <- function(Y, X, ncomp){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  PCs <- pca$x[, 1:ncomp]
  PCs <- as.matrix(PCs)
  
  fit <- stats::.lm.fit(x = PCs, y = Y)
  Y.hat <- PCs%*%fit$coefficients
  
  return(Y.hat)
}


