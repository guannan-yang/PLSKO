# source('/data/gpfs/projects/punim0613/Guan/Sim2/RSource/plsko_adv.R')
# source('/data/gpfs/projects/punim0613/Guan/Sim2/RSource/plsko_plsonly.R')
# source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/spca.R")
# source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/internal.R")
# source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/knockoff_filters.R")
# source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/simknockoffs.R")
# source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/KOBT/create.pc.knockoff.R")

# source('./RSource/plsko_adv.R')
source('C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/plsko_plsonly.R')
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/spca.R")
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/seqknockoff/internal.R")
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/seqknockoff/knockoff_filters.R")
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/seqknockoff/simknockoffs.R")
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/KOBT/create.pc.knockoff.R")
source("C:/Users/guannany/OneDrive - The University of Melbourne/PhD/vitual_desk/Sim2/RSource/KO_perf.R")

# source('./RSource/plsko_plsonly.R')
# source("./RSource/spca.R")
# source("./RSource/seqknockoff/internal.R")
# source("./RSource/seqknockoff/knockoff_filters.R")
# source("./RSource/seqknockoff/simknockoffs.R")
# source("./RSource/KOBT/create.pc.knockoff.R")


# Part II 1) fdr: calculating fdr

fdr <- function(S, beta.true) {
  fdp = sum(beta.true[S] == 0)/max(1, length(S))
  return(fdp)
}

fdr.factor <- function(S, true.sig){
  fdp = sum(!(S %in% true.sig))/max(1, length(S))
  return(fdp)
}
###########################################################################

###########################################################################
# Part II: additional functions
# Part II 2) pow: calculating power

pow <- function(S, beta.true) {
  tdp = sum(beta.true[S] != 0)/sum(beta.true != 0)
  return(tdp)
}

pow.factor <- function(S, true.sig){
  tdp = sum(true.sig %in% S)/max(1, length(true.sig))
  return(tdp)
}

###########################################################################

###########################################################################
# n: number of observations
# p: number of variables

# DGPX.AR1 <- function(n,p,r=5,theta=1, method = "AR1", rho =0.5, X.dis = "Normal"){
#   if (method =="AR1" & X.dis == "Normal") {
#     obj <- norm.sim(n, p, rho = rho)
#     
#   } else if (method == "half" & X.dis=="Normal") {
#     obj <- norm.half.sim(n, p, rho = rho)
#     
#   } else if (method == "AR1" & X.dis == "NB"){
#     sim <- neg.bin.sim(n, p, rho = rho)
#     X <- sim$X
#     # mimic preprocess
#     ## adjust library size
#     X <- X/matrix(rep(rowMeans(X),p),n, byrow = F)
#     ## log2 transform
#     X <- log2(X+0.1)
#     
#     obj <- list(X = X, cov = sim$cov)
#   } else if (method == "half" & X.dis == "NB"){
#     
#     sim <- neg.bin.half.sim(n, p, rho = rho)
#     X <- sim$X
#     # mimic preprocess
#     ## adjust library size
#     X <- X/matrix(rep(rowMeans(X),p),n, byrow = F)
#     ## log2 transform
#     X <- log2(X+0.1)
#     
#     obj <- list(X = X, cov = sim$cov)
#   }
#   return(obj)
# }
# 
# # copula test: can it simulate negative binominal
# 
# neg.bin.sim <- function(n, p, rho, c.mu = 400, c.phi = 0.1){
#   Sigma <- diag(1,nrow = p, ncol = p)
#   for (i in 1:p) {
#     for (j in 1:p) {
#       Sigma[i,j] <- rho^(abs(i-j))
#     }
#   }
#   mu <- rep(0, p)
#   X <- MASS::mvrnorm(n, mu=mu, Sigma=Sigma)
#   u <- pnorm(X)
#   
#   size = 1/c.phi
#   #qnbinom(p, size, prob, mu, lower.tail = TRUE, log.p = FALSE)
#   X.negbin <- qnbinom(u, size=size, mu = c.mu)
#   X.negbin[is.na(X.negbin)] <- 1
#   return(list(X = X.negbin,
#               cov = Sigma))
# }

norm.sim <- function(n, p, rho){
  Sigma <- diag(1,nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  mu <- rep(0, p)
  X <- MASS::mvrnorm(n, mu=mu, Sigma=Sigma)
  u <- pnorm(X)
  
  X.norm <- qnorm(u)
  return(list(X = X.norm,
              cov = Sigma))
}

norm.half.sim <- function(n, p, rho){
  Sigma <- diag(1,nrow = p, ncol = p)
  for (i in 1:round(p/2)) {
    for (j in 1:round(p/2)) {
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  
  mu <- rep(0, p)
  X <- MASS::mvrnorm(n, mu=mu, Sigma=Sigma)
  u <- pnorm(X)
  
  X.norm <- qnorm(u)
  return(list(X = X.norm,
              cov = Sigma))
}

neg.bin.half.sim <- function(n, p, rho, c.mu = 400, c.phi = 0.1){
  Sigma <- diag(1,nrow = p, ncol = p)
  for (i in 1:round(p/2)) {
    for (j in 1:round(p/2)) {
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  mu <- rep(0, p)
  X <- MASS::mvrnorm(n, mu=mu, Sigma=Sigma)
  u <- pnorm(X)
  
  size = 1/c.phi
  #qnbinom(p, size, prob, mu, lower.tail = TRUE, log.p = FALSE)
  X.negbin <- qnbinom(u, size=size, mu = c.mu)
  X.negbin[is.na(X.negbin)] <- 1
  return(list(X = X.negbin,
              cov = Sigma))
}
###########################################################################

###########################################################################
# n: number of observations
# p: number of variables
# s: level of sparsity
# A: amplitude of signals
# method: 1- linear: a 1 x n vector y is generated form normal distribution
#         2- nonlinear:  a 1 x n vector y is generated form logistic distribution

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

###########################################################################
# Knockoff function
###########################################################################

knockoffX.sim <- function(X, method, r = 5, ..., ko.diagnosis = F){
  # source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/plsko_adv.R")
  # source('/data/gpfs/projects/punim0613/Guan/Sim2/RSource/plsko_plsonly.R')
  # source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/internal.R")
  # source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/knockoff_filters.R")
  # source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/seqknockoff/simknockoffs.R")
  # source('./RSource/hdko.R')
  
  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)
  
  ncomp = r
  
  start.time <- Sys.time()
  
  if(method == "default"){
    cov <- cov(X)
    ko <- create.gaussian(X, mu, cov)
  }
  else if(method == "JS"){
    cov <- matrix(as.numeric(corpcor::cov.shrink(X,verbose=F)), nrow=ncol(X))
    ko <- create.second_order(X, shrink = T)
  }
  else if(method == "glasso"){
    X.glasso <- CVglasso::CVglasso(X = X, diagonal = F)
    glasso.cov <- pre2cov.glasso.chol(X, X.glasso)
    
    cov <- glasso.cov
    ko <- create.gaussian(X, mu, cov, method = "sdp")
  }
  else if(method == "lorec"){
    # tuning?
    lorec <- lorec(crossprod(X), lambda=.01, delta=0.01)
    lorec.cov <- lorec$L + lorec$S
    
    cov <- lorec.cov
    ko <- create.gaussian(X, mu, cov)
  }

  else if(method == "poet"){
    cov <- poet.cov(X, ncomp = r)
    ko <- create.gaussian(X, mu, cov)
  }
  else if(method == "ipad.pca.center"){
    cov = NA
    ko <- IPAD_ko(X, r = r)
  }

  else if(method == "KOScreen"){
    cov <- NA
    ko <- create.MK(X, 1:p, M = 1, method = 'svd.full', bigmemory = F)[[1]]
  }
  else if(method == "PLSKO"){
    cov <- NA
    ko <- plsko.adv(X, ncomp = r, ...)
    ko <- ko$X_k
  }
  else if(method == "PLSKO.only"){
    cov <- NA
    ko <- plsko.nools(X, ncomp = r, ...)
    ko <- ko$X_k
  }
  else if(method == "knn.ko"){
    cov <- NA
    ko <- knn.ko(X)
  }
  
  else if(method == "seqknockoffs"){
    cov <- NA
    X <- as.data.frame(X)
    X[is.na(X)] <- 0
    ko <- knockoffs_seq(X)
  }
  
  else if(method == "perm"){
    cov <- NA
    ko <- X[sample(n),]
  }
  else if(method == "PCKO"){
    cov <- NA
    ko <- hd.ko(X, fit.method = "pca", r.fit = r, ...)
  }
  
  end.time <- Sys.time()
  run.time <- difftime(end.time, start.time, units = "secs")
  
  if(ko.diagnosis){
  diagnosis <- dis.diagnosis(X, ko, mrc.ncomp = r) 
  result <- list(#cov = cov,
                  ko = ko,
                  method = method,
                  run.time = run.time,
                  diagnosis = diagnosis
  )
  }
  else{result <- list(#cov = cov,
                      ko = ko,
                      method = method,
                      run.time = run.time
  )}
  
  return(result)
}


####################
# Knockoff filter
####################
ko.filter <- function(X, Xk, y, q, method = "lasso.lcd", sparsity = 0.1, data.log = T, ...){
  n = nrow(X)
  p = ncol(X)

  if(method == "lasso.lcd"){
    W <- stat.lasso_coefdiff(X, Xk, y, ...)
  }
  else if(method == "lasso.max.lambda"){
    W <- stat.lasso_lambdadiff(X, Xk, y, ...)
  }
  else if(method == "pls.regression.lcd"){
    W <- W_spls_lcd(X, Xk, y, lcd_comp = 1)
  }
  
  else if(method == "plsda.lcd"){ # exclude when continuous y design
    W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1)
  }
  else if(method == "spls.regression.lcd"){
    keepX <- round(ncol(X)*sparsity)+1
    W <- W_spls_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  }
  else if(method == "splsda.lcd"){ # exclude when continuous y design
    keepX <- round(ncol(X)*sparsity)+1
    W <- W_splsda_lcd(X, Xk, y, ncomp = 2, lcd_comp = 1, keepX = rep(keepX, 2))
  }
  
  else if(method == "lasso.logistic"){ # exclude when continuous y design
    W <- stat.lasso_coefdiff_bin(X, Xk, y, ...)
  }
  else if(method == "DEG"){
    W <- W_DEG(X, Xk, y, data.log = data.log)
  }
  else if(method == "RF"){
    W <- stat.random_forest(X, Xk, y)
  }
  
  # find the knockoff threshold T
  t = sort(c(0, abs(W)))
  
  ratio = c(rep(0, p))
  ratio.plus = c(rep(0, p))
  
    for (j in 1:p) {
      ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
      ratio.plus[j] = (1+sum(W <= -t[j]))/max(1, sum(W >= t[j]))
    }

  id = which(ratio <= q)[1]
  id.plus = which(ratio.plus <= q)[1]
  
  if(length(id) == 0){
    T = Inf
  } else {
    T = t[id]
  }
  
  if(length(id.plus) == 0){
    T.plus = Inf
  } else {
    T.plus = t[id.plus]
  }
  # set of discovered variables
  S = integer0_test(which(W >= T))
  S.plus = integer0_test(which(W >= T.plus))
  
  return(list(S = S, S.plus = S.plus, W = W))
}

ko.filter.mv <- function(X, Xk, y, q, method, rank = 5, plus = F, adjust.q = F, sparse = 0.1, add = F, refer.group =1){
  n = nrow(X)
  p = ncol(X)
  
  if(method == "plsda.lcd"){
    W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1:rank)
  }
  if(method == "splsda.lcd"){
    keepX = rep(round(sparse * p)+1, rank)
    W <- W_splsda_lcd(X, Xk, y, lcd_comp = 1:rank, keepX = keepX)
  }
  if(method == "lasso.mn"){
    W <- W_mnlasso_lcd(X, Xk, y, refer.group = refer.group)
  }
  
  S = list()
  
  rank = ncol(W)
  for(i.comp in 1:rank){
    # find the knockoff threshold T for each class/component
    t = sort(c(0, abs(W[,i.comp])))
    ratio = c(rep(0, p))
    
    if(!plus){
      for (j in 1:p) {
        ratio[j] = (sum(W[,i.comp] <= -t[j]))/max(1, sum(W[,i.comp] >= t[j]))
      }
    }
    else{
      for (j in 1:p) {
        ratio[j] = (1+sum(W[,i.comp] <= -t[j]))/max(1, sum(W[,i.comp] >= t[j]))
      }
    }
    id = which(ratio <= q)[1]
    if(length(id) == 0){
      T = Inf
    } else {
      T = t[id]
    }
    
    # set of discovered variables
    S[[i.comp]] = integer0_test(which(W[,i.comp] >= T))
  }
  
  return(S)
}

############################################################
############################################################ 

W_splsda_lcd <- function(X, X_k, y, ncomp = 5, lcd_comp = 1:ncomp, keepX = rep(ncol(X),ncomp)){
  #colnames(X) <- as.character(1:ncol(X))
  #colnames(X_k) <- paste0(as.character(1:ncol(X_k)),"k")
  nvar <- ncol(X)
  
  #swap X and X_k column randomly
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  # X.swap  = X * (1-swap.M) + X_k * swap.M
  # Xk.swap = X * swap.M + X_k * (1-swap.M)
  X_X_k.swap <- cbind(X*(1-swap.M)+ X_k*swap.M, X*swap.M + X_k*(1-swap.M))
  colnames(X_X_k.swap) <- as.character(1:ncol(X_X_k.swap))
  
  splsda <- mixOmics::splsda(X_X_k.swap, y, ncomp = ncomp, keepX = keepX)
  #plotIndiv(plsda, legend = T)
  
  splsda_lcd <- abs(splsda$loadings$X[1:nvar, lcd_comp])- abs(splsda$loadings$X[(nvar+1):(nvar*2), lcd_comp])
  splsda_lcd <- unname(splsda_lcd)
  
  #swap back
  swap.lcd <- matrix(swap, nrow = ncol(X), ncol = length(lcd_comp), byrow = F)
  splsda_lcd <- splsda_lcd * (1-2*swap.lcd)
  
  return(splsda_lcd)}



W_spls_lcd <- function(X, X_k, y, ncomp = 5, lcd_comp = 1:ncomp, keepX = rep(ncol(X),ncomp)){
  #colnames(X) <- as.character(1:ncol(X))
  #colnames(X_k) <- paste0(as.character(1:ncol(X_k)),"k")
  nvar <- ncol(X)
  
  #swap X and X_k column randomly
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  # X.swap  = X * (1-swap.M) + X_k * swap.M
  # Xk.swap = X * swap.M + X_k * (1-swap.M)
  X_X_k.swap <- cbind(X*(1-swap.M)+ X_k*swap.M, X*swap.M + X_k*(1-swap.M))
  colnames(X_X_k.swap) <- as.character(1:ncol(X_X_k.swap))
  
  spls <- mixOmics::spls(X_X_k.swap, y, ncomp = ncomp, keepX = keepX)
  #plotIndiv(plsda, legend = T)
  
  spls_lcd <- abs(spls$loadings$X[1:nvar, lcd_comp])- abs(spls$loadings$X[(nvar+1):(nvar*2), lcd_comp])
  spls_lcd <- unname(spls_lcd)
  
  #swap back
  swap.lcd <- matrix(swap, nrow = ncol(X), ncol = ncomp, byrow = F)
  spls_lcd <- spls_lcd * (1-2*swap.lcd)
  
  return(spls_lcd)}

W_DEG <- function(X, Xk, y, data.log = T){
  X <- t(X)
  Xk <- t(Xk)
  if(data.log ) X.bind <- rbind(2^X, 2^Xk)
  else {X.bind <- rbind(X, Xk)
  X.bind[X.bind<0] <- 0}
  DGE.list <- DGEList(counts = X.bind, group = as.factor(y))
  design <- model.matrix(~DGE.list$samples$group)
  v <- voom(DGE.list, design)
  fit <- lmFit(v, design)
  efit <- eBayes(fit)
  t <- efit$t
  W <- abs(t)[1:nrow(X),2] - abs(t)[(nrow(X)+1):(nrow(X)*2),2]
  return(W)
}

W_mnlasso_lcd <- function(X, X_k, y, refer.group =1){
  
  nvar <- ncol(X)
  
  #swap X and X_k column randomly
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  # X.swap  = X * (1-swap.M) + X_k * swap.M
  # Xk.swap = X * swap.M + X_k * (1-swap.M)
  X_X_k.swap <- cbind(X*(1-swap.M)+ X_k*swap.M, X*swap.M + X_k*(1-swap.M))
  colnames(X_X_k.swap) <- as.character(1:ncol(X_X_k.swap))
  
  coefs <- cv_coeffs_glmnet(X_X_k.swap, y, family = "multinomial", parallel = F)
  Z = data.frame(row.names = as.character(1:(2*nvar)))
  for (i.class in 1:length(coefs)){
    Z[,i.class] <- abs(coefs[[refer.group]][2:(2*nvar+1)]) + abs(coefs[[i.class]][2:(2*nvar+1)])
  } 
  # remove reference group's coef
  Z <- Z[,-refer.group]
  
  W = abs(Z[1:nvar,]) - abs(Z[(nvar+1):(2*nvar),])
  
  swap.lcd <- matrix(swap, nrow = ncol(X), ncol = length(coefs)-1, byrow = F)
  
  W = W * (1-2*swap.lcd)
  return(W)
}

#### 
cv_coeffs_glmnet <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  # Standardize variables
  X = scale(X)

  n = nrow(X); p = ncol(X)

  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family

  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, intercept=intercept,
                                     standardize=F,standardize.response=F, parallel=parallel, ...)

  coef(cv.glmnet.fit, s = "lambda.min")
}

integer0_test <- function(which){
  if(identical(which, integer(0))){
    return(NULL)
  }
  else {return(which)}
}

fFact <- function(X, r){
  obj <- RSpectra::svds(X,r)
  if(r>1){
  C <- obj$u %*% diag(obj$d) %*% t(obj$v)  
  }
  else{
    C <- obj$u %*% t(obj$v) *obj$d
  }
  return(C)
}
IPAD_ko <- function(X, r = 5){
  n = nrow(X)
  p = ncol(X)
  C = fFact(X,r)
  res = X - C
  resScale = sqrt(apply(res^2,2,sum))
  resScale = diag(resScale)
  e = matrix(rnorm(n*p),n,p)
  eScale = sqrt(apply(e^2,2,sum))
  eScale = diag(eScale^-1)
  nres = e %*% resScale %*% eScale
  X_k = C + nres
  return(X_k)
}

pre2cov.glasso <- function(glasso){
  pre <- glasso$Omega
  cov <- solve(pre)
  cov.sym <- (cov+t(cov))/2
  return(cov.sym)
}

pre2cov.glasso.chol <- function(X,glasso){
  pre <- glasso$Omega
  chol <- chol(pre)
  inv.pre <- chol2inv(chol)
  # scaling
  emp.cov <- cov(X)
  r <- sqrt(diag(emp.cov)/diag(inv.pre))
  cov <- diag(r)%*%inv.pre%*%diag(r)
  
  return(cov)
}

S_BH <- function(X, y, q = 0.05){
  raw.p <- apply(X,2, function(Xj){
    p <- summary(lm(Xj~y))[["coefficients"]][2,4]
  })
  
  adj.p <- p.adjust(raw.p, method = "BH")
  S <- which(adj.p < q)
  return(S)
}

jaccard <- function(a, b){
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


r_criterion <- function(X, rmax = 20){
  # rmax = 10 # maximum number of factors
  n = nrow(X)
  p = ncol(X)
  
  PC = matrix(NA,rmax+1,1)
  SX = scale(X)
  mXX = SX %*%t(SX)
  for(k in rmax:0){
    if (k == 0){
      PC[k+1,] = sum(SX^2/(n*p))
    } else 
      eig <- RSpectra::eigs_sym(mXX, k, which = "LM", sigma = NULL, lower = TRUE)
    meigvec <- eig$vectors 
    mF <- sqrt(n) * meigvec       # estimated factors
    Lam <- (t(mF) %*% SX)/n
    if(k==rmax){
      sigma2 = sum((SX-mF %*% Lam)^2)/(n*p)
    }
    PC[k+1,] = sum((SX-mF %*% Lam)^2)/(n*p) + k*sigma2*((n+p)/(n*p))*log((n*p)/(n+p))
  }
  r = which(PC == min(PC))[1]-1
}

hdknockoff <- function(X,y, method = "PLSKO.only", ..., screen.method = "margin", screen.family = "gaussian", p.thres = 0.1, keep.num = 3*ncol(X)){
  n = nrow(X)
  p = ncol(X)
  sample0 <- sample(round(n/2))
  sample1 <- base::setdiff(1:n, sample0)
  # y0 = y[1:round(n/2)]
  # X0 = X[1:round(n/2),]
  # y1 = y[(round(n/2)+1):n]
  # X1 = X[(round(n/2)+1):n,]
  y0 = y[sample0]
  X0 = X[sample0,]
  y1 = y[sample1]
  X1 = X[sample1,]
  
  if(screen.method == "EN"){
    Op_lambda <- cv.glmnet(X0,y0, family = screen.family, alpha = 0.5)$lambda.min
    # apply glmnet to find beta vector under original scale
    Beta <- glmnet(X0,y0,lambda = Op_lambda, family = screen.family, alpha = 0.5)$beta
    Beta <- as.vector(Beta)
    screened <- Beta!=0
  }
  else if(screen.method == "lasso"){
    #Finding Optimal value of Lambda
    Op_lambda <- cv.glmnet(X0,y0, family = screen.family)$lambda.min
    # apply glmnet to find beta vector under original scale
    Beta <- glmnet(X0,y0,lambda = Op_lambda, family = screen.family)$beta
    Beta <- as.vector(Beta)
    screened <- Beta!=0
  }
  else if(screen.method=="margin"){
    #use marginal regression to screen variable: raw.p < 0.1
    raw.p <- apply(X0,2, function(Xj){
      p <- summary(lm(Xj~y0))[["coefficients"]][2,4]
    })
    screened <- (raw.p < p.thres) & (rank(raw.p) <= keep.num)
  }

  r_hat <- r_criterion(X1[, screened], rmax = 10)
  X_ko <- knockoffX.sim(X1[, screened],method = method, r = r_hat, ..., ko.diagnosis = F)$ko
  X_ko <- rbind(X0[, screened], X_ko)
  X_ko <- X_ko[rownames(X),] #recover sample order
  Xnew = X[, screened]
  return(list(Xnew = Xnew,
              X_ko = X_ko))
}
