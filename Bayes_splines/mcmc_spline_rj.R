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
  
  p_beta <- function(sig_sq, tau_sq = 1000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    #browser()
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    rmvnorm(1, mean = mu, sigma = sig)
  }
  
  bhat <- function(sig_sq, tau_sq = 10000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    mu
  }
  
  met_gibbs <- function(its, max = 50) {
    
    ### RJ MCMC Ratio Calculations ###
    
    post <- function(tau_sq = 0.01, g1 = 10001, g2 = 10000, Xcurr = X_curr, Xcand = X_cand, .Y = y, n = length(y)) {
      ts <- 1 / tau_sq
      V_curr <- solve(t(Xcurr) %*% Xcurr + ts * diag(ncol(Xcurr)))
      V_cand <- solve(t(Xcand) %*% Xcand + ts * diag(ncol(Xcand)))
      
      num <- 0.5 * log(det(ts * diag(ncol(V_curr)))) + 0.5 * log(det(V_cand))
      den <- 0.5 * log(det(ts * diag(ncol(V_cand)))) + 0.5 * log(det(V_curr))
      
      ahatcurr <- bhat(sig_sq = 1, X = Xcurr)
      ahatcand <- bhat(sig_sq = 1, X = Xcand)
      d_curr <- g2 + (t(.Y) %*% .Y) - (t(ahatcurr) %*% solve(V_curr) %*% ahatcurr)
      d_cand <- g2 + (t(.Y) %*% .Y) - (t(ahatcand) %*% solve(V_cand) %*% ahatcand)
      d_term <- (g1 + (n/2)) * (log(d_curr/d_cand))
      result <- (num - den) + d_term
      if(is.na(result)) browser()
      return(result) #simplifying everything to logs
    }
    
    p <- function(lambda = 1, .k) {
      num <- .k * log(lambda)
      den <- log(exp(lambda) - 1) + lfactorial(.k)
      return(num - den)
    }
    
    prior <- function(k) {
      return(0.5 * (p(.k = k+1) - p(.k = k)))
    }
    
    proposal <- function(k, c = 0.4) {
      d <- 1/3#c * min(0, log(p(.k = k+1) - p(.k = k)))
      b <- 1/3#c * min(0, log(p(.k = k) - p(.k = k+1)))
      num <- log(d) - log(k)
      den <- log(b) - log(2)
      return(num - den)
    }
    
    ## Initializing Values for RJMCMC ##
    
    nknot <- 0 # First iteration must be a birth
    X_cand <- matrix(NA, length(x), max)
    X_cand[,1] <- 1
    mat_t <- matrix(NA, its, max)
    mat_s <- matrix(NA, its, max)
    mat_t <- rep(0, max)
    mat_s <- rep(1, max)
    
    X_curr <- X_cand

    mat_beta <- matrix(NA, its, max)
    mat_beta <- rep(1, ncol(X_curr))
    mat_sig <- rep(NA, its)
    mat_sig <- 0.1
    
    
    for(it in 2:its) {
      
      # need to generate number of knots and produce X_cand
      # you must choose birth in first iteration
      # sample knot, sign
      # if accept, you have 1 BF
      # next step: B, D, or C
      
      propose_t <- function(x) {
        log(dunif(x, 0, 1))
      } #for dealing with proposal knot locations
      
      samp <- function(knots = nknot) {
        if(knots == 0) {return(1)}
        if(knots == max) {sample(2:3, 1)}
        sample(3, 1)
        # 1 = BIRTH
        # 2 = DEATH
        # 3 = CHANGE
      }
      choice <- samp()
      
      if(choice == 1) {
        basis_vec <- rep(NA, length(x))
        candidate_t <- runif(1)
        candidate_s <- sample(c(-1,1), 1, replace = TRUE)
        
        pos <- function(v) {
          v[v < 0] <- 0
          v
        }
        
        basis_vec <- pos(candidate_s * (x - candidate_t))
        
        X_cand <- cbind(X_curr, basis_vec)
        
        Xcand_revised <- X_cand[,colSums(is.na(X_cand)) != nrow(X_cand)] %>% 
          as.matrix()
        Xcurr_revised <- X_curr[,colSums(is.na(X_curr)) != nrow(X_curr)] %>% 
          as.matrix()
        # if a column is all NAs, then that number must equal the number of rows

        ratio_rj <- post(Xcurr = Xcurr_revised, Xcand = Xcand_revised) + prior(k = nknot+1) + proposal(k = nknot+1)
        
        accept_prob <- min(0, ratio_rj)
        if(is.na(ratio_rj)) browser()
        if(log(runif(1)) < accept_prob) { 
          nknot <- nknot + 1
          mat_t <- rbind(mat_t, candidate_t)
          mat_s <- rbind(mat_s, candidate_s)
          X_curr <- X_cand
        } else {
          mat_t[it,] <- mat_t[it-1,]
          mat_s[it,] <- mat_s[it-1,]
        }
      } 
      
      else if(choice == 2) { #DEATH
        pick <- sample(2:(nknot+1), 1)
        X_curr <- X_curr[,-pick]
        nknot <- nknot - 1
      } 
      
      else { #CHANGE
        X_cand <- X_curr
        col <- sample(2:(nknot+1), 1)
        #row <- sample(1:length(x), 1)
        basis_vec <- rep(NA, length(x))
        candidate_t <- runif(1)
        candidate_s <- sample(c(-1,1), 1, replace = TRUE)
        pos <- function(v) {
          v[v < 0] <- 0
          v
        }
        basis_vec <- pos(candidate_s * (x - candidate_t))
        X_cand[,col] <- basis_vec  
        
        ratio_change <- function(tau_sq = 0.01, g1 = 0.01, g2 = 0.01, n = length(y)) {
          Vprime <- solve(t(X_cand) %*% X_cand + tau_sq * diag(ncol(X_cand)))
          V <- solve(t(X_curr) %*% X_curr + tau_sq * diag(ncol(X_curr)))
          aprime <- Vprime %*% t(X_cand) %*% y 
          a <- V %*% t(X_curr) %*% y 
          dprime <- g2 + (t(y) %*% y) - (t(aprime) %*% solve(Vprime) %*% aprime)
          d <- g2 + (t(y) %*% y) - (t(a) %*% solve(V) %*% a)
          V_part <- 0.5 * log(det(Vprime)) -  0.5 * log(det(V))
          d_part <- (g1 + n/2) * (log(d) - log(dprime))
          V_part + d_part
        }
        
        acc <- min(0, ratio_change())
        
        if(log(runif(1)) < acc) {
          X_curr <- X_cand
        }
      }
      
      X_curr2 <- X_curr[,colSums(is.na(X_curr)) != nrow(X_curr)]
      #browser()
      mat_beta <- rbind(mat_beta, p_beta(sig_sq = mat_sig[it-1], X = X_curr2))
      #browser()
      mat_sig <- c(mat_sig, p_sig(X = X_curr2, beta = mat_beta[it,]))
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
