logsumchoose <- function(q, vec) {
  #Calculate log choose via log gamma functions
  logchoose <- lgamma(q+1) - lgamma(q-vec+1) - lgamma(vec+1)
  
  #Find the maximum element among the log choose's
  maxlogchoose <- max(logchoose)
  
  #Numeric trick for calculating log sum choose
  est <- maxlogchoose + log(sum(exp(logchoose - maxlogchoose)))
  return(est)
}

SoftThreshold <- function(x, lambda) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function(sigma, i, mu, 
                                maxiter = 50, threshold = 1e-2) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol = 1.5, 
                          mu = NULL, maxiter = 50, threshold = 1e-2, 
                          verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    print(i)
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

remMap.grid <- function(X, Y, standardize = TRUE, 
                        n.lambda = 10, lambda.min.ratio = 0.01) {
  
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  #Find maximum inner product of X_i and Y_j
  lambda.vec <- rep(NA, q)
  for (j in 1:q) {
    lambda.vec[j] <- max(abs(colSums(X*Y[,j])))
  }
  lambda.max <- max(lambda.vec)
  
  #Compute lambda sequence
  lambda.seq <- round(exp(seq(log(lambda.max), log(lambda.max*lambda.min.ratio), 
                              length.out = n.lambda)), digits = 10)
  remMap.lambda.grid <- expand.grid(lambda.seq, lambda.seq)
  colnames(remMap.lambda.grid) <- c("lambda1", "lambda2")
  
  return(remMap.lambda.grid)
}

remMap.one <- function(X, Y, standardize = TRUE, 
                       lambda1, lambda2, C = NULL) {
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  if (is.null(C)) {
    Theta0 <- t(remMap(X, Y, lambda1, lambda2)$phi)
  } else {
    Theta0 <- t(remMap(X, Y, lambda1, lambda2, C.m = t(C))$phi)
  }
  
  return(Theta0)
}

remMap.EBIC <- function(X, Y, Theta0, standardize = TRUE, 
                        gamma = 1) {
  
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  E <- Y - X %*% t(Theta0)
  RSS <- diag(t(E) %*% E)
  EBIC <- n * sum(log(RSS)) + 
    log(n) * sum(Theta0 != 0) + 
    2 * gamma * logsumchoose(q, rowSums(t(Theta0) != 0))
  
  return(EBIC)
}

#remMap all in one function
remMap.whole <- function(X, Y, standardize = TRUE, 
                         n.lambda = 10, lambda.min.ratio = 0.01, 
                         C = NULL, gamma = 1) {
  
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  #Find maximum inner product of X_i and Y_j
  lambda.vec <- rep(NA, q)
  for (j in 1:q) {
    lambda.vec[j] <- max(abs(colSums(X*Y[,j])))
  }
  lambda.max <- max(lambda.vec)
  
  #Compute lambda sequence
  lambda.seq <- round(exp(seq(log(lambda.max), log(lambda.max*lambda.min.ratio), 
                              length.out = n.lambda)), digits = 10)
  remMap.lambda.grid <- expand.grid(lambda.seq, lambda.seq)
  colnames(remMap.lambda.grid) <- c("lambda1", "lambda2")
  
  n.lambda.sq <- dim(remMap.lambda.grid)[1]
  
  #List for storing Theta0 matrices for each (lambda1, lambda2) pair
  ls <- vector("list", n.lambda.sq)
  
  for (i in 1:n.lambda.sq) {
    print(i)
    if (is.null(C)) {
      ls[[i]] <- t(remMap(X, Y, remMap.lambda.grid[i,1], remMap.lambda.grid[i,2])$phi)
    } else {
      ls[[i]] <- t(remMap(X, Y, remMap.lambda.grid[i,1], remMap.lambda.grid[i,2], C.m = t(C))$phi)
    }
  }
  
  #Calculate EBIC for all grids
  EBIC.vec <- rep(NA, n.lambda.sq)
  
  for (i in 1:n.lambda.sq) {
    E <- Y - X %*% t(ls[[i]])
    RSS <- diag(t(E) %*% E)
    EBIC.vec[i] <- n * sum(log(RSS)) + 
      log(n) * sum(ls[[i]] != 0) + 
      2 * gamma * logsumchoose(q, rowSums(t(ls[[i]]) != 0))
  }
  
  #Grid that minimizes EBIC
  opt.idx <- which.min(EBIC.vec)
  
  #Return Theta0 chosen using EBIC and corresponding (lambda1, lambda2)
  return(list(Theta0 = ls[[opt.idx]], 
              lambda1.opt = remMap.lambda.grid[opt.idx,1],
              lambda2.opt = remMap.lambda.grid[opt.idx,2]))
}

precM.glasso <- function(X, standardize = TRUE, 
                         gamma = 0.5, rholist = NULL,
                         thr = 1e-4, maxit = 1e4) {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Standardize X to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
  }
  
  CovM <- t(X) %*% (X) / n
  
  #Array of precision matrix candidates
  a <- glassopath(CovM, rholist, thr, maxit)
  
  #Calculate EBIC for all grids
  n.grid <- dim(a$wi)[3]
  EBIC.vec <- rep(NA, n.grid)
  
  for (i in 1:n.grid) {
    L <- n * (tr(CovM %*% a$wi[,,i]) - determinant(a$wi[,,i])$modulus[1]) / 2
    E <- sum(a$wi[,,i] != 0)
    EBIC.vec[i] <- 2 * L + E * log(n) + 4 * gamma * E * log(p)
  }
  
  #Grid that minimizes EBIC
  opt.idx <- which.min(EBIC.vec)
  
  #Return precM chosen using EBIC
  return(precM = a$wi[,,opt.idx])
}

#Should standardize X to avoid calling standardization for p times
precM.NL.one <- function(X, row.idx, method = "EBIC", gamma = 0.5,
                         n.lambda = 100, lambda = NULL) {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #use cross-validation for selecting row
  if (method == "CV") {
    cv.fit <- cv.glmnet(x = X[,-row.idx], y = X[,row.idx], 
                        family = "gaussian", nlambda = n.lambda,
                        lambda = lambda, standardize = FALSE, 
                        intercept = FALSE)
    lambda.opt <- cv.fit$lambda[cv.fit$index[2]] #Lambda selected
    beta <- cv.fit$glmnet.fit$beta[,cv.fit$index[2]] #Corresponding lasso coefficients
    SSE <- colSums((X[,row.idx,drop=FALSE] - X[,-row.idx] %*% beta)^2)
    tau.sq <- SSE / n + lambda.opt * sum(abs(beta))
  } else {
    #Use EBIC for selecting row (default)
    fit <- glmnet(x = X[,-row.idx], y = X[,row.idx], 
                  family = "gaussian", nlambda = n.lambda,
                  lambda = lambda, standardize = FALSE, 
                  intercept = FALSE)
    lambda <- fit$lambda #Lambda grid used in fitting lasso
    n.grid <- length(lambda) #Number of grid points
    SSE.vec <- colSums((matrix(X[,row.idx], n, n.grid) 
                        - X[,-row.idx] %*% fit$beta)^2)
    df.vec <- colSums(fit$beta != 0) #Number of non-zero parameters in coefficients across grid
    EBIC.vec <- SSE.vec + df.vec * log(n) + 2 * gamma * lchoose(p, df.vec)
    opt.idx <- which.min(EBIC.vec)
    beta <- fit$beta[,opt.idx]
    tau.sq <- SSE.vec[opt.idx] / n + lambda[opt.idx] * sum(abs(fit$beta[,opt.idx]))
  }
  ls <- list(beta = beta, tau.sq = tau.sq)
  return(ls)
}

precM.NL.whole <- function(X, standardize = TRUE, method = "EBIC",
                           gamma = 0.5, n.lambda = 100) {
  
  p <- dim(X)[2]
  
  #Standardize X to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
  }
  
  C.hat <- diag(1, p, p)
  tau.sq.vec <- rep(NA, p)
  
  for (i in 1:p) {
    ls <- precM.NL.one(X, i, method, gamma, n.lambda)
    C.hat[i,-i] <- -ls$beta
    tau.sq.vec[i] <- ls$tau.sq
  }
  
  precM <- diag(1/tau.sq.vec) %*% C.hat
  return(precM)
}

precM.QO <- function(X, standardize = TRUE) {
  
  n <- dim(X)[1]
  
  #Standardize X to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
  }
  
  CovM <- t(X) %*% (X) / n
  
  precM <- InverseLinfty(CovM, n, resol = 1.3)
  return(precM)
}

precM <- function(X, method = "glasso", standardize = TRUE) {
  
  #Nodewise lasso
  if (method == "NL") {
    precM <- precM.NL.whole(X, standardize)
  } else 
    #Quadratic optimization
    if (method == "QO") {
      precM <- precM.QO(X, standardize)
    } else {
      #Use glasso as default
      precM <- precM.glasso(X, standardize)
    }
  return(precM)
}

#Calculate grid for DrFARM for a given number of latent factor k
DrFARM.grid <- function(X, Y, Theta0, precM, k,
                        lambda1.opt, lambda2.opt, 
                        standardize = TRUE) {
  
  n <- dim(X)[1]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  E.star <- Y - X %*% Theta0.db.t
  fa.res <- fa(E.star, nfactors = k, rotate = "none", scores = "regression", fm = "ml", covar = TRUE)
  diag.Psi <- fa.res$uniquenesses
  
  Psi.range <-  unname(summary(diag.Psi)[-4])
  DrFARM.lambda.grid <- expand.grid(lambda1.opt / Psi.range, lambda2.opt / Psi.range)
  
  colnames(DrFARM.lambda.grid) <- c("lambda1", "lambda2")
  
  return(DrFARM.lambda.grid)
  
}

DrFARM.one <- function(X, Y, Theta0, precM, k, 
                       lambda1, lambda2, C = NULL,
                       standardize = TRUE,
                       thres = 1e-4,
                       rotate = "none",
                       scores = "regression",
                       fm = "ml") {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  #Initialization
  
  #Debias the sparse estimator used for initial value
  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  #Calculate the residual matrix E.star
  E.star <- Y - X %*% Theta0.db.t
  
  #Obtain the factor analysis component using fa() from psych
  fa.res <- fa(E.star, nfactors = k, rotate = rotate, scores = scores, fm = fm, covar = TRUE)
  B <- fa.res$loadings
  diag.Psi <- fa.res$uniquenesses
  Psi <- diag(diag.Psi)
  Psi.Inv <- diag(1/diag.Psi)
  PsiInv.B <- Psi.Inv %*% B
  Bt.PsiInv.B <- t(B) %*% Psi.Inv %*% B
  InvList <- list()
  #To be converted into Rcpp?
  for (i in 1:n) {
    sandInv <- solve(diag(1/d[i], k, k) + Bt.PsiInv.B)
    InvList[[i]] <- Psi.Inv - PsiInv.B %*% sandInv %*% t(PsiInv.B)
  }
  
  iter <- 0
  prev.loss <- 1e300
  diff <- 1e300
  
  Theta.t <- t(Theta0)
  prev.Theta.t <- Theta.t
  prev.B <- B
  prev.Psi <- diag.Psi
  
  start <- Sys.time()
  while (diff > thres) {
    iter <- iter + 1
    e <- Y - X %*% Theta.t
    #E-step
    eet <- t(e) %*% e
    #sum of all expected zzt's
    # E.zzt <- n * (diag(k) - W %*% B) + W %*% eet %*% t(W)
    # E.zzt.inv <- solve(E.zzt)
    
    WList <- list()
    for (i in 1:n) {
      WList[[i]] <- d[i] * t(B) %*% InvList[[i]]
    }
    #sum of all expected z's
    E.zt <- matrix(NA, k, n)
    for (i in 1:n) {
      E.zt[,i] <- WList[[i]] %*% e[i,]
    }
    #checkZ[[iter]] <- E.zt
    
    E.zzt.List <- list()
    for (i in 1:n) {
      E.zzt.List[[i]] <- d[i] * (diag(1, k, k) - WList[[i]] %*% B) + E.zt[,i] %*% t(E.zt[,i])
    }
    
    E.zzt <- Reduce("+", E.zzt.List)
    E.zzt.inv <- solve(E.zzt)
    
    #i <- as.numeric(task_id)
    Y.aug <- Y - t(B %*% E.zt)
    if (is.null(C)) {
      Theta.t <- remMap(X.m = X, Y.m = Y.aug, lamL1 = lambda1, lamL2 = lambda2, phi0 = prev.Theta.t, C.m = NULL, sigma = diag.Psi)$phi
    } else {
      Theta.t <- remMap(X.m = X, Y.m = Y.aug, lamL1 = lambda1, lamL2 = lambda2, phi0 = prev.Theta.t, C.m = t(C), sigma = diag.Psi)$phi
    }
    e <- Y - X %*% Theta.t
    Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n
    E.star <- Y - X %*% Theta.db.t #Bia-corrected residual matrix
    
    B <- t(E.zt %*% E.star) %*% E.zzt.inv
    EET.star <- t(E.star) %*% E.star
    
    diag.Psi <- diag(EET.star - B %*% E.zt %*% E.star)/n
    Psi.Inv <- diag(1/diag.Psi)
    PsiInv.B <- Psi.Inv %*% B
    Bt.PsiInv.B <- t(B) %*% Psi.Inv %*% B
    
    InvList <- list()
    for (i in 1:n) {
      sandInv <- solve(diag(1/d[i], k, k) + Bt.PsiInv.B)
      InvList[[i]] <- Psi.Inv - PsiInv.B %*% sandInv %*% t(PsiInv.B)
    }
    #i <- as.numeric(task_id)
    loss <- sum(e^2 * matrix(1/diag.Psi, n, q, byrow = TRUE)) / 2 +
      lambda1 * sum(abs(Theta.t)) +
      lambda2 * sum(sqrt(rowSums(Theta.t^2))) +
      n * sum(log(diag.Psi)) / 2
    #n * determinant(B %*% t(B) + Psi)$modulus[1] / 2
    print(c("M-step", "iter", iter, loss))
    #print(loss)
    diff <- prev.loss - loss
    prev.loss <- loss
    if (diff > 0) {
      prev.Theta.t <- Theta.t
      prev.B <- B
      prev.Psi <- diag.Psi
    } else if (diff < 0) {
      iter <- iter - 1
      break
    }
    #}
    Theta.t <- prev.Theta.t
    B <- prev.B
    diag.Psi <- prev.Psi
    
    #checkloss[[iter]] <- loss
    prev.loss <- loss
  }
  return(list(Theta = t(Theta.t), B = B, E.Z = t(E.zt), diag.Psi = diag.Psi))
}

DrFARM.EBIC <- function(X, Y, Theta, B, E.Z, diag.Psi, standardize = TRUE, 
                        gamma = 1) {
  
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  E <- Y - X %*% t(Theta) - E.Z %*% t(B)
  
  ratio <- diag(t(E) %*% E) / diag.Psi
  EBIC <- sum(ratio) + n * sum(log(diag.Psi)) +
    log(n) * sum(Theta != 0) + 
    2 * gamma * logsumchoose(q, rowSums(t(Theta) != 0))
  
  return(EBIC)
}

DrFARM.whole <- function(X, Y, Theta0, precM, k, 
                         lambda1.opt, lambda2.opt, C = NULL,
                         standardize = TRUE,
                         gamma = 1,
                         thres = 1e-4,
                         rotate = "none",
                         scores = "regression",
                         fm = "ml") {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  Theta0.db.t <- t(Theta0) + t(t(Y - X %*% t(Theta0)) %*% X %*% precM) / n
  E.star <- Y - X %*% Theta0.db.t
  fa.res <- fa(E.star, nfactors = k, rotate = "none", scores = "regression", fm = "ml", covar = TRUE)
  diag.Psi <- fa.res$uniquenesses
  
  Psi.range <-  unname(summary(diag.Psi)[-4])
  DrFARM.lambda.grid <- expand.grid(lambda1.opt / Psi.range, lambda2.opt / Psi.range)
  
  colnames(DrFARM.lambda.grid) <- c("lambda1", "lambda2")
  
  n.lambda.sq <- dim(DrFARM.lambda.grid)[1]
  
  #List for storing Theta0 matrices for each (lambda1, lambda2) pair
  ls <- vector("list", n.lambda.sq)
  
  for (i in 1:n.lambda.sq) {
    print(i)
    ls[[i]] <- DrFARM.one(X, Y, Theta0, precM, k, 
                          DrFARM.lambda.grid[i,1], DrFARM.lambda.grid[i,2],
                          C = C, standardize = FALSE,
                          thres = thres, rotate = rotate,
                          scores = scores, fm = fm)
  }
  
  #Calculate EBIC for all grids
  EBIC.vec <- rep(NA, n.lambda.sq)
  
  for (i in 1:n.lambda.sq) {
    E <- Y - X %*% t(ls[[i]]$Theta) - ls[[i]]$E.Z %*% t(ls[[i]]$B)
    
    ratio <- diag(t(E) %*% E) / ls[[i]]$diag.Psi
    EBIC.vec[i] <- sum(ratio) + n * sum(log(ls[[i]]$diag.Psi)) +
      log(n) * sum(ls[[i]]$Theta != 0) + 
      2 * gamma * logsumchoose(q, rowSums(t(ls[[i]]$Theta) != 0))
  }
  
  #Grid that minimizes EBIC
  opt.idx <- which.min(EBIC.vec)
  
  #Return Theta0 chosen using EBIC and corresponding (lambda1, lambda2)
  return(list(Theta = ls[[opt.idx]]$Theta,
              B = ls[[opt.idx]]$B,
              E.Z = ls[[opt.idx]]$E.Z,
              diag.Psi = ls[[opt.idx]]$diag.Psi,
              lambda1.opt = DrFARM.lambda.grid[opt.idx,1],
              lambda2.opt = DrFARM.lambda.grid[opt.idx,2]))
}

entry.pvalue <- function(X, Y, Theta, B, E.Z, precM, 
                         standardize = TRUE) {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  Y.aug <- Y - E.Z %*% t(B)
  
  Theta.t <- t(Theta)
  Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n
  
  s <- colSums(Theta.t != 0) #degree of freedom by trait
  sses <- diag(t(Y - X %*% Theta.t - E.Z %*% t(B))  %*% (Y - X %*% Theta.t - E.Z %*% t(B)))
  Psi.star <- sses / (n - s)
  
  CovM <- t(X) %*% X / n
  sqrtPhi <- sqrt(diag(precM %*% CovM %*% t(precM)))
  
  Z <- matrix(NA, p, q)
  for (i in 1:q) {
    Z[,i] <- sqrt(n) * Theta.db.t[,i] / (sqrt(Psi.star[i]) * sqrtPhi)
  }
  
  pval <- 2 * pnorm(-abs(Z))
  return(pval)
}

pleio.pvalue <- function(X, Y, Theta, B, E.Z, precM, 
                         standardize = TRUE) {
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  #Standardize both X and Y to mean 0 and variance 1
  if (standardize == TRUE) {
    X <- scale(X)
    Y <- scale(Y)
  }
  
  Y.aug <- Y - E.Z %*% t(B)
  
  Theta.t <- t(Theta)
  Theta.db.t <- Theta.t + t(t(Y.aug - X %*% Theta.t) %*% X %*% precM) / n
  
  s <- colSums(Theta.t != 0) #degree of freedom by trait
  sses <- diag(t(Y - X %*% Theta.t - E.Z %*% t(B))  %*% (Y - X %*% Theta.t - E.Z %*% t(B)))
  Psi.star <- sses / (n - s)
  
  CovM <- t(X) %*% X / n
  sqrtPhi <- sqrt(diag(precM %*% CovM %*% t(precM)))
  
  Z <- matrix(NA, p, q)
  for (i in 1:q) {
    Z[,i] <- sqrt(n) * Theta.db.t[,i] / (sqrt(Psi.star[i]) * sqrtPhi)
  }
  
  pval <- 2 * pnorm(-abs(Z))
  
  T <- rowSums(tan(pi * (0.5 - pval))) / q
  pval2 <- 2 * pcauchy(-abs(T))
  return(pval2)
}
