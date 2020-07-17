library(mvtnorm)
mcmc_spline <- function(nknot = 5) {
  
  p_t <- function(sig_sq, betahat, X) {
    (-1/(2*sig_sq)) * t(y - (X %*% betahat)) %*% (y - (X %*% betahat))
  }
  
  p_sig <- function(a = 1, b = 1, X, beta) {
    n <- nrow(X)
    a_term <- a + (n/2)
    b_term <- 0.5 * (2*b + (t(y - (X %*% beta)) %*% (y - (X %*% beta))))
    1 / rgamma(1, shape = a_term, rate = b_term)
  }
  
  p_beta <- function(sig_sq, tau_sq = 10000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    rmvnorm(1, mean = mu, sigma = sig)
  }
  
  bhat <- function(sig_sq, tau_sq = 10000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    mu
  }
  
  met_gibbs <- function(its) {
    
    mat_t <- matrix(NA, its, nknot)
    mat_t[1,] <- seq(0.1, 0.5, by = 0.1)
    
    ## INITIALIZE X_curr ##
    X_curr <- matrix(NA, nknot, length(x))
    s <- rep(1, nknot)
    for(i in 1:nknot) {
      for(j in 1:length(x)) {
        X_curr[i,j] <- max(s[i] * (x[j] - mat_t[1,i]), 0)
      } 
    }
    X_curr <- t(X_curr)
    X_curr <- cbind(1, X_curr)
    
    mat_beta <- matrix(NA, its, ncol(X_curr))
    mat_sig <- rep(NA, its)
    
    mat_sig[1] <- 0.1
    mat_beta[1,] <- rep(1, ncol(X_curr))
    
    ar <- 0
    
    for(it in 2:its) {
      propose_t <- function(x) {
        log(dunif(x, 0, 1))
      }
      
      candidate_t <- runif(nknot, 0, 1)
      
      ### CREATE BASIS FUNCTIONS FROM THIS t-vector  ###
      
      s <- rep(1, nknot)
      
      X_cand <- t(X_curr[,-1]) #matrix(NA, nknot, length(x))
      
      for(i in 1:nknot) {
        for(j in 1:length(x)) {
          X_cand[i,j] <- max(s[i] * (x[j] - candidate_t[i]), 0)
        } #creating basis
      } #X Matrix
      
      X_cand <- t(X_cand)
      X_cand <- cbind(1, X_cand)
      
      ratio <- 
        p_t(sig_sq = mat_sig[it-1], betahat = bhat(mat_sig[it-1], X = X_cand), X = X_cand) -
        p_t(sig_sq = mat_sig[it-1], betahat = bhat(mat_sig[it-1], X = X_curr), X = X_curr)
      
      accept_prob <- min(0, ratio)
      u <- log(runif(1))
      #browser()
      if(u < accept_prob) {
        ar <- ar + 1
        mat_t[it,] <- candidate_t
        X_curr <- X_cand
      } else {
        mat_t[it,] <- mat_t[it-1,]
      }
      
      mat_beta[it,] <- p_beta(sig_sq = mat_sig[it-1], X = X_curr)
      mat_sig[it] <- p_sig(beta = mat_beta[it,], X = X_curr)
    }
    list(mat_beta, mat_sig, mat_t, X_curr, ar/its)
  }
  a <- met_gibbs(its = its)
  mat_beta <- a[[1]] #beta values
  mat_sig <- a[[2]] #sig^2 values
  mat_t <- a[[3]]
  X_curr <- a[[4]]
  ar <- a[[5]]
  
  list(mat_beta, mat_sig, mat_t, X_curr, ar)
}


### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(nknot = 5, knots) {
  s = rep(1, nknot)
  Xm <- matrix(NA, nknot, length(x))
  for(i in 1:nknot) {
    for(j in 1:length(x)) {
      Xm[i,j] <- max(s[i] * (x[j] - knots[i]), 0)
    } 
  } 
  Xm <- t(Xm)
  Xm <- cbind(1, Xm)
  Xm
}
