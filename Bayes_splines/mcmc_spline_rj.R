library(mvtnorm)
library(dplyr)
mcmc_spline <- function() {
  
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
    
    ### RJ MCMC Ratio Calculations ###
    
    post <- function(tau_sq = 0.01, g1 = 0.01, g2 = 0.01, Xcurr = X_curr, Xcand = X_cand, .Y = y) {
      ts <- 1 / tau_sq
      V_curr <- solve(t(Xcurr) %*% Xcurr + ts * diag(ncol(Xcurr)))
      V_cand <- solve(t(Xcand) %*% Xcand + ts * diag(ncol(Xcand)))
      
      num <- 0.5 * log(det(ts * diag(ncol(V_curr)))) + 0.5 * log(det(V_cand))
      dem <- 0.5 * log(det(ts * diag(ncol(V_cand)))) + 0.5 * log(det(V_curr))
      
      ahatcurr <- bhat(sig_sq = 1, X = Xcurr)
      ahatcand <- bhat(sig_sq = 1, X = Xcand)
      d_curr <- g2 + (t(.Y) %*% .Y) - (t(ahatcurr) %*% solve(V_curr) %*% ahatcurr)
      d_cand <- g2 + (t(.Y) %*% .Y) - (t(ahatcand) %*% solve(V_cand) %*% ahatcand)
      d_term <- g1 * (log(d_curr) - log(d_cand))
      
      return((num - den) + d_term) #simplifying everything to logs
    }
    
    p <- function(lambda = 1, .k) {
      num <- lambda ^ .k
      den <- (exp(lambda) - 1) * factorial(.k)
      return(num/den)
    }
    
    prior <- function(k) {
      return(0.5 * (p(.k = k+1) - p(.k = k)))
    }
    
    proposal <- function(k, c = 0.4) {
      d <- c * min(0, log(p(.k = k+1) - p(.k = k)))
      b <- c * min(0, log(p(.k = k) - p(.k = k+1)))
      num <- log(d) - log(k)
      den <- log(b) - log(2)
      return(num - den)
    }
    
    ## First MCMC Iteration ##
    
    nknot <- 1 # First iteration must be a birth
    X_cand <- matrix(NA, nknot, length(x))
    
    for(i in 1:nknot) {
      for(j in 1:length(x)) {
        X_cand[i,j] <- max(s[i] * (x[j] - mat_t[1,][i]), 0)
      } 
    }
    
    X_cand <- t(X_cand)
    X_cand <- cbind(1, X_cand)
    X_curr <- X_cand
    # The first step (birth) done manually
    
    mat_t <- matrix(NA, its, nknot)
    mat_s <- matrix(NA, its, nknot)
    mat_t[1,] <- rep(0, nknot)
    mat_s[1,] <- rep(1, nknot)
    
    mat_beta <- matrix(NA, its, ncol(X_curr))
    mat_sig <- rep(NA, its)
    mat_sig[1] <- 0.1
    mat_beta[1,] <- rep(1, ncol(X_curr))
    
    for(it in 2:its) {
      
      # need to generate number of knots and produce X_cand
      # you must choose birth in first iteration
      # sample knot, sign
      # if accept, you have 1 BF
      # next step: B, D, or C
      
      propose_t <- function(x) {
        log(dunif(x, 0, 1))
      } #for dealing with proposal knot locations
      
      samp <- function() {
        sample(3, 1)
        # 1 = BIRTH
        # 2 = DEATH
        # 3 = CHANGE
      }
      choice <- samp()
      
      if(choice == 1) {
        nknot <- nknot + 1
        basis_vec <- rep(NA, nknot)
        candidate_t <- runif(nknot, 0, 1)
        candidate_s <- sample(c(-1,1), nknot, replace = TRUE)
        
        for(j in seq_len(length(x))) {
          basis_vec[j] <- max(candidate_s * (x[j] - candidate_t), 0)
        }
        X_cand <- cbind(X_cand, basis_vec)
        
        ratio_rj <- post() + prior(k = nknot) + proposal(k = nknot)
        
        accept_prob <- min(0, ratio_rj)
        
        if(log(runif(1)) < accept_prob) {
          mat_t[it,] <- candidate_t
          mat_s[it,] <- candidate_s
          X_curr <- X_cand
        } else {
          mat_t[it,] <- mat_t[it-1,]
          mat_s[it,] <- mat_s[it-1,]        
        }
      } 
      
      else if(choice == 2) {
        pick <- sample(2:nknot, 1)
        X_cand <- X_cand[,-pick]
        X_curr <- X_cand
        nknot <- nknot - 1
      } 
      
      else {
        col <- sample(2:nknot, 1)
        row <- sample(1:length(x), 1)
        #X_cand[row, col] <- choose a new basis function value based on a new knot location
      }
      
      mat_beta[it,] <- p_beta(sig_sq = mat_sig[it-1], X = X_curr)
      mat_sig[it] <- p_sig(beta = mat_beta[it,], X = X_curr)
    }
    list(mat_beta, mat_sig, mat_t, X_curr, ar/its, mat_s)
  }
  a <- met_gibbs(its = its)
  mat_beta <- a[[1]] #beta values
  mat_sig <- a[[2]] #sig^2 values
  mat_t <- a[[3]]
  X_curr <- a[[4]]
  ar <- a[[5]]
  mat_s <- a[[6]]
  
  list(mat_beta, mat_sig, mat_t, X_curr, ar, mat_s)
}


### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(nknot, knots, signs) {
  
  Xm <- matrix(NA, nknot, length(x))
  for(i in 1:nknot) {
    for(j in 1:length(x)) {
      Xm[i,j] <- max(signs[i] * (x[j] - knots[i]), 0)
    } 
  } 
  Xm <- t(Xm)
  Xm <- cbind(1, Xm)
  Xm
}
