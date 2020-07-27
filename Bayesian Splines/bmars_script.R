library(mvtnorm)
library(dplyr)

pos <- function(vec) {
  ((abs(vec) + vec) / 2)
} # set anything less than 0 to 0; creating basis functions [a]+

mcmc_spline <- function(its, max_knot = 50, max_j = 3) {
  
  p_t <- function(sig_sq, betahat, X) {
    (-1/(2*sig_sq)) * t(y - (X %*% betahat)) %*% (y - (X %*% betahat))
  } # full conditional for the knot locations, t
  
  p_sig <- function(a = 1, b = 1, X, beta) {
    n <- nrow(X)
    a_term <- a + (n/2)
    b_term <- 0.5 * (2*b + (t(y - (X %*% beta)) %*% (y - (X %*% beta))))
    1 / rgamma(1, shape = a_term, rate = b_term)
  } # full conditional for sigma^2 in the regression equation
  
  p_beta <- function(sig_sq, tau_sq = 1000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    rmvnorm(1, mean = mu, sigma = sig)
  } # full conditional for the regression coefficients
  
  bhat <- function(sig_sq, tau_sq = 10000, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau_sq) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    mu
  } #marginalized betahat for regression coefficients
  # THIS IS THE GAUSSIAN LIKELIHOOD
  
  met_gibbs <- function(its) {
    
    ### RJ MCMC Conditionals ###
    
    post <- function(tau_sq = 0.0001, g1 = 10001, g2 = 10000, Xcurr = X_curr, Xcand = X_cand, .Y = y) {
      n <- length(y)
      ts <- 1 / tau_sq
      V_curr <- solve(t(Xcurr) %*% Xcurr + tau_sq * diag(ncol(Xcurr)))
      V_cand <- solve(t(Xcand) %*% Xcand + tau_sq * diag(ncol(Xcand)))
      
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
    } #RJMCMC posterior
    
    p <- function(lambda = 1, .k) {
      num <- .k * log(lambda)
      den <- log(exp(lambda) - 1) + lfactorial(.k)
      return(num - den)
    } #RJMCMC poisson prior with default lambda = 1
    
    prior_b <- function(k) {
      return(0.5 * (p(.k = k+1) - p(.k = k)))
      # needs to be (0.5) ^ J_k+1
    } #prior in RJMCMC ratio for a birth step
    
    prior_d <- function(k) {
      return(0.5 * (p(.k = k-1) - p(.k = k)))
      # needs to be (0.5) ^ J_k+1
    } #prior in RJMCMC ratio for a death step
    
    proposal_b <- function(k) {
      d <- 1/3 #probability of death
      b <- 1/3 #probability of birth
      # NOTE: P(change) = 1/3
      num <- log(d) - log(k)
      den <- log(b) - log(2) 
      # should be log(b) - J_k+1 * log(2*n)
      return(num - den)
    } #proposal in RJMCMC ratio for a birth step
    
    proposal_d <- function(k) {
      d <- 1/3 #probability of death
      b <- 1/3 #probability of birth
      # NOTE: P(change) = 1/3
      num <- log(d) - log(k)
      den <- log(b) - log(2)
      return(den - num)
    } #proposal in RJMCMC ratio for a death step
    
    ## Initializing Values for RJMCMC ##
    
    nknot <- 0 # First iteration must be a birth
    X_curr <- rep(1, length(x)) %>% as.matrix() 
    # first current X-matrix is just an intercept
    
    mat_t <- array(NA, dim = c(its, max_knot, max_j))
    mat_s <- array(NA, dim = c(its, max_knot, max_j))
    mat_v <- array(NA, dim = c(its, max_knot, max_j))
    mat_j <- matrix(NA, its, max_knot)
    
    mat_beta <- matrix(NA, its, max_knot)
    mat_beta[1,1] <- mean(y)
    mat_sig <- rep(NA, its)
    mat_sig[1] <- 1

    
    for(it in 2:its) {
      
      fate <- function(knots = nknot) {
        if((knots == 0) | (knots == 1)) {return(1)} # having 0 or 1 knots must auto defer to birth
        if(knots == max_knot) {sample(2:3, 1)} #at max_knot knot capacity, can only delete or change
        sample(3, 1)
        # 1 = BIRTH
        # 2 = DEATH
        # 3 = CHANGE
      }
      
      samp_j <- function(limit = max_j) {
        sample(1:limit, 1)
      }
      
      choice <- fate()
      
      if(choice == 1) { # BIRTH
        # sample j
        # sample j random knots, j random signs, j variables
        # make 1 basis function from this
        # ALMOST the same
        # Prior and proposal will be different
        # evaluate likelihood
        # accept/reject
        j <- samp_j()
        
        samp_x <- function(vars = j) {
          n <- ncol(X)
          sample(2:n, j) %>% sort() #accounting for intercept term, sort the terms
        } #ARE THESE THE v_ij TERMS
        
        cols <- samp_x() # NEED TO ACCOUNT FOR THE TYPE
        candidate_t <- runif(j, 0, 1) #change this
        candidate_s <- sample(c(-1,1), j, replace = TRUE)
        Xmat <- X[,cols]
        
        unsign <- t(t(Xmat) - candidate_t)
        basis_mat <- t(candidate_s * t(unsign)) 
        basis_mat <- pos(basis_mat) # matrix of basis functions, same size as X_mat
        
        # MULTIPLY BASIS FUNCTIONS TOGETHER TO GET THE DIFFERENT "TYPES"
        
        X_cand <- cbind(X_curr, basis_mat)
        
        ratio_rj <- post(Xcurr = X_curr, Xcand = X_cand) + prior_b(k = nknot) + proposal_b(k = nknot)
        #like * prior / proposal
        accept_prob <- min(0, ratio_rj)
        
        if(is.na(ratio_rj)) browser()
        
        if(log(runif(1)) < accept_prob) { 
          nknot <- nknot + 1
          mat_t[it,] <- mat_t[it-1,]
          mat_t[it,nknot] <- candidate_t
          mat_s[it,] <- mat_s[it-1,]
          mat_s[it,nknot] <- candidate_s
          X_curr <- X_cand
        } else {
          mat_t[it,] <- mat_t[it-1,]
          mat_s[it,] <- mat_s[it-1,]
        }
      } 
      
      else if(choice == 2) { #DEATH
        # select function to delete
        # prior, proposal different
        pick <- sample(2:(nknot+1), 1)
        X_cand <- X_curr[,-pick]
        
        ratio_rj <- post(Xcurr = X_curr, Xcand = X_cand) + prior_d(k = nknot) + proposal_d(k = nknot)
        
        accept_prob <- min(0, ratio_rj)
        if(is.na(ratio_rj)) browser()
        
        if(log(runif(1)) < accept_prob) { 
          nknot <- nknot - 1
          X_curr <- X_cand
          mat_t[it,(1:nknot)] <- mat_t[it-1,(1:(nknot+1))[-(pick-1)]]
          mat_s[it,(1:nknot)] <- mat_s[it-1,(1:(nknot+1))[-(pick-1)]]
        } else{
          mat_t[it,] <- mat_t[it-1,]
          mat_s[it,] <- mat_s[it-1,]         
        }
      } 
      
      else { #CHANGE
        X_cand <- X_curr
        col <- sample(2:(nknot+1), 1)
        #row <- sample(1:length(x), 1)
        candidate_t <- runif(1)
        candidate_s <- sample(c(-1,1), 1, replace = TRUE)
        
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
          mat_t[it,] <- mat_t[it-1,]
          mat_t[it, (col-1)] <- candidate_t
          mat_s[it,] <- mat_s[it-1,]
          mat_s[it, (col-1)] <- candidate_s
        } else{
          mat_t[it,] <- mat_t[it-1,]
          mat_s[it,] <- mat_s[it-1,]         
        }
      }

      mat_beta[it,1:(nknot+1)] <- p_beta(sig_sq = mat_sig[it-1], X = X_curr)
      mat_sig[it] <- p_sig(X = X_curr, beta = mat_beta[it,(1:(nknot+1))])
    }
    list(mat_beta, mat_sig, mat_t, mat_s, X_curr)
  }
  a <- met_gibbs(its = its)
  mat_beta <- a[[1]] #beta values
  mat_sig <- a[[2]] #sig^2 values
  mat_t <- a[[3]]
  mat_s <- a[[4]]
  X_curr <- a[[5]]
  
  list(mat_beta, mat_sig, mat_t, mat_s, X_curr)
}


### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(nknot, knots, signs) {
  
  Xm <- matrix(NA, nknot, length(x))
  for(i in 1:nknot) {
    for(j in 1:length(x)) {
      Xm[i,j] <- pos(signs[i] * (x[j] - knots[i]))
    } 
  } 
  Xm <- t(Xm)
  Xm <- cbind(1, Xm)
  Xm
}
