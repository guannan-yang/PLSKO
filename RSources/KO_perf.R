gram.mat <- function(X, X_k){
  
  # scale
  X <- scale(X)
  X_k <- scale(X_k)
  
  gram_xx <- cor(X)
  gram_xkxk <- cor(X_k, X_k)
  gram_xxk <- cor(X, X_k)
  
  s_hat <- diag(gram_xxk)
  # generate a matrix with col as c(i,j, Gxx, Gxkxk, Gxxk, i==j)
  gram_xx.melt <- gram_xx %>% melt() %>% rename(gram_xx = value)
  gram_xkxk.melt <- gram_xkxk %>% melt()%>% rename(gram_xkxk = value)
  gram_xxk.melt <- gram_xxk %>% melt()%>% rename(gram_xxk = value)
  
  gram.melt <- gram_xx.melt %>% 
    cbind(gram_xkxk = gram_xkxk.melt$gram_xkxk) %>% 
    cbind(gram_xxk = gram_xxk.melt$gram_xxk) %>% 
    mutate(ij = (Var1==Var2))
  
  return(gram.melt)
}

draw.gram.mat <- function(gram.melt, method1.name = "", method2.name = ""){
  
  # too many points, randomly select a part of rows to draw plots
  if(nrow(gram.melt)> 1e+05) sampled.row <- sample.int(nrow(gram.melt), 1000)
  else sampled.row <- 1:nrow(gram.melt)
  
  gram.melt.sub <- rbind(gram.melt[sampled.row,], gram.melt[gram.melt$ij==T,])
  
  scatter1 <- ggplot(gram.melt.sub, aes(x = gram_xx, y = gram_xkxk))+
    geom_point(aes(color = ij))+
    xlab("cor(X_i,X_j)")+
    ylab("cor(Xk_i,Xk_j)")+
    ggtitle(paste0(method1.name, " " ,method2.name))+
    labs(color = "Corresponding knockoff (i==j)")
  print(scatter1)
  
  scatter2 <- ggplot(gram.melt.sub, aes(x = gram_xx, y = gram_xxk))+
    geom_point(aes(color = ij))+
    xlab("cor(X_i,X_j)")+
    ylab("cor(X_i,Xk_j)")+
    ggtitle(paste0(method1.name, " ", method2.name))+
    labs(color = "Corresponding knockoff (i==j)")
  print(scatter2)
}

cov.diagno <- function(Z1, Z2){
  tr1.tr2 <- 0
  tr1.2 <- 0
  n <- nrow(Z1)
  for(i in 1:n){
    for(j in c(1:n)){
      if(i!=j){tr1.tr2 <- tr1.tr2 + crossprod(Z1[i,],Z1[j,])^2+ crossprod(Z2[i,],Z2[j,])^2}
      else{tr1.2 <- tr1.2 + crossprod(Z1[i,],Z2[j,])^2}
    }
  }
  phi.cov <- 1/(n*(n-1))*tr1.tr2 - 2/n^2*tr1.2
  phi.cov <- phi.cov[1,1]
  phi.cov
}

# MMD by kernlab


# KNN1
knn1.diagno <- function(Z1,Z2){
  rownames(Z1) <- paste0(rownames(Z1),"_1")
  rownames(Z2) <- paste0(rownames(Z2),"_2")
  n <- nrow(Z1)
  
  Z <- rbind(Z1,Z2)
  dist.matrix <- as.matrix(dist(Z))
  # exclude the self distance and the distance with its own knockoffs
  for(i in 1:n){
    dist.matrix[i,i] <- NA
    dist.matrix[i, i+n] <- NA
  }
  for(i in (n+1):(2*n)){
    dist.matrix[i,i] <- NA
    dist.matrix[i, i-n] <- NA
  }
  
  nn <- apply(dist.matrix,2,function(x)return(rownames(dist.matrix)[which.min(x)]))
  
  nn.i1 <- (names(nn) %in% rownames(Z1))&(nn %in% rownames(Z1))
  nn.i2 <- (names(nn) %in% rownames(Z2))&(nn %in% rownames(Z2))
  
  phi.nn <- 1/(2*n)*(sum(nn.i1) + sum(nn.i2))
  return(phi.nn)
}

# Energy
energy.phi <- function(Z1,Z2){
  n <- nrow(Z1)
  norm12 <- 0
  norm11 <- 0
  norm22 <- 0
  normii <- 0
  
  for (i in 1:n) {
    for (j in 1:n) {
      norm12 <- norm12 + sqrt(crossprod(Z1[i,]-Z2[j,]))
      norm11 <- norm11 + sqrt(crossprod(Z1[i,]-Z1[j,]))
      norm22 <- norm22 + sqrt(crossprod(Z2[i,]-Z2[j,]))
    }
    normii <- normii + sqrt(crossprod(Z1[i,]-Z2[i,]))
  }
  
  dist.energy <- 2/n^2*norm12 - 1/n^2*(norm11+norm22)
  phi.energy <- n/2*(dist.energy - 2/n^2*normii)
  phi.energy <- phi.energy[1,1]
  
  
  return(phi.energy)
}

# second-order gram based distance 
so.gram <- function(X, Xk){
  p <- ncol(X)
  X <- scale(X, scale = F)
  X_k <- scale(X, scale = F)
  gram.xx <- cor(X)
  gram.xkxk <- cor(Xk)
  gram.xxk <- cor(X, Xk)
  M = matrix(1, p, p)
  diag(M) <- 0
  so = norm(gram.xx-gram.xxk, type = "F") +
    norm(M * (gram.xx - gram.xxk), type = "F")
  
  decorr <- norm(diag(gram.xxk), type = '2')
  
  return(list(so = so, decorr = decorr))
}


swap.cbind <- function(X, X_k){
  # avoiding duplicate names
  # rownames(X_k) <- paste0(rownames(X_k),"_k")
  # colnames(X_k) <- paste0(colnames(X_k),"_k")
  
  # center 
  X <- apply(X, 2, function(x){x-mean(x)})
  X_k <- apply(X, 2, function(x){x-mean(x)})
  XXk <- cbind(X,X_k)
  
  # swapping half of columns of X and X_k
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X_k * (1-swap.M)+ X * swap.M
  XXk.swap <- cbind(X.swap, Xk.swap)
  
  # avoid duplicate names
  # rownames(XXk.swap) <- paste0(rownames(XXk.swap),"_swap")
  # colnames(XXk.swap) <- paste0(colnames(XXk.swap),"_swap")  
  
  return(list(XXk = XXk,
              XXk.swap = XXk.swap))
}

L.mrc <- function(X, Xk, res.method = "linear", ncomp = NULL){
  p = ncol(X)
  n = nrow(X)
  colnames(X) <- 1:p
  colnames(Xk) <- paste0(1:p,"k")
  
  l.mrc <- 0
  
  for (j in 1:p){
    X.j <- X[,j]
    X._j <- X[,-j]
    
    if(res.method =="linear"){
      fit <- linear.regression.generator(Y = X.j ,X = cbind(X._j, Xk))
    }
    else if(res.method =="pls"){
      X.run = cbind(X._j, Xk)
      
      if(is.null(ncomp)){
        ncomp <- ceiling(min(n/2, p))+1
      }
      
      keepX <-  rep(2*p-2, ncomp)
      fit <- plsda.recovery.generator(X.j, X.run, ncomp = ncomp, keepX = keepX)
    }
    else if(res.method == "block.pls"){
      # warning(colnames(X._j))
      # warning(colnames(Xk))

      X.run <- list(X._j = X._j, Xk = Xk)
      Y <- cbind(X.j, X.j)
      Y <- as.matrix(Y)
      colnames(Y) <- c(1,2)
      rownames(Y) <- rownames(X)
      
      if(is.null(ncomp)){
        ncomp <- ceiling(min(n/2, p/2))+1 #make sure this.ncomp > 2
      }

      keepX <- list(X._j = rep(p-1, ncomp),
                    Xk = rep(p, ncomp))
      # full design
      des.mat <- matrix(1, nrow = 3, ncol = 3)
      diag(des.mat) <- 0
      colnames(des.mat) <- c("X._j", "Xk", "Y")
      rownames(des.mat) <- c("X._j", "Xk", "Y")
      
      fit <- block.pls.recovery.generator(Y, X.run, ncomp = ncomp, design = des.mat)
    }
    
    res <- X.j - fit
    res.var <- var(res)
    l.mrc <- l.mrc + 1/res.var
  }
  
  return(l.mrc)
}

dis.diagnosis <- function(X, Xk, mrc.ncomp = round(min(dim(X))/2)){
  
  Xswap <- swap.cbind(X, Xk)
  XXk <- Xswap$XXk
  XXk.swap <- Xswap$XXk.swap
  # covariance diagnosis
  dis.cov <- cov.diagno(XXk,XXk.swap)
  
  # MMD (Radial Basis Gaussian kernel)
  dis.mmd <- max(kmmd(XXk, XXk.swap, kernel="rbfdot")@mmdstats)
  
  # knn1 diagnosis
  # dis.knn1 <- knn1.diagno(XXk,XXk.swap)
  
  # energy diagonosis
  dis.energy <- energy.phi(XXk, XXk.swap)
  
  
  # gram matrix (w/o swap)
  corr.gram <- so.gram(X, Xk)$so
  decorr.gram <- so.gram(X, Xk)$decorr 
  
  #l.mrc <- L.mrc(X, Xk, res.method = "pls", ncomp = mrc.ncomp)
  
  return(list(dis.cov = dis.cov,
              dis.mmd = dis.mmd,
              # dis.knn1 = dis.knn1,
              dis.energy = dis.energy,
              gram.corr = corr.gram,
              gram.decorr = decorr.gram#,l.mrc = l.mrc
              ))
}

