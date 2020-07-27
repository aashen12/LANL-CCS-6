library(mvtnorm)
library(dplyr)

### [a]+ FUNCTION ### 

pos <- function(vec) {
  ((abs(vec) + vec) / 2)
} # set anything less than 0 to 0; creating basis functions [a]+

#######################################################################################

### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(nknot, knots, signs, data, deg = 1) {
  tdat <- t(data)
  temp1 <- pos(signs * (tdat[vars,,drop = F] - knots))^degree
  
  if(length(vars) == 1) {
    return(c(temp1))
  } else {
    temp2 <- 1
    for(pp in 1:length(vars)) {
      temp2 <- temp2 * temp1[pp,]
    }
    return(temp2)
  }
}

#######################################################################################

bmars <- function(its=3000, max_knot = 50, max_j = 3, tau_sq = 10^4, g1=0, g2=0, h1=10, h2=10) {
  Xt <- t(X)
  n <- length(y)
  p <- ncol(X)
  ssy <- sum(y^2)
  p_t <- function(sig_sq, betahat, X) {
    (-1/(2*sig_sq)) * t(y - (X %*% betahat)) %*% (y - (X %*% betahat))
  } # full conditional for the knot locations, t
  
  p_sig <- function(a = 1, b = 1, X, beta) {
    n <- nrow(X)
    a_term <- a + (n/2)
    b_term <- 0.5 * (2*b + (t(y - (X %*% beta)) %*% (y - (X %*% beta))))
    1 / rgamma(1, shape = a_term, rate = b_term)
  } # full conditional for sigma^2 in the regression equation
  
  p_beta <- function(sig_sq, tau_sq = 10000, X) {
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
    
  ## Initializing Values for RJMCMC ##

  
  mat_t <- array(NA, dim = c(its, max_knot, max_j)) #knots
  mat_s <- array(NA, dim = c(its, max_knot, max_j)) #signs
  mat_v <- array(NA, dim = c(its, max_knot, max_j)) #vars
  mat_j <- matrix(NA, its, max_knot) #number of interactions
  
  mat_beta <- matrix(NA, its, max_knot + 1)
  mat_beta[1,1] <- mean(y)
  mat_sig <- lam <- rep(NA, its)
  mat_sig[1] <- 1
  lam[1] <- 1
  
  nknot <- 0 # First iteration must be a birth
  X_curr <- rep(1, n) %>% as.matrix() 
  # first current X-matrix is just an intercept

  Vinv_curr <- crossprod(X_curr) + 1/tau_sq
  ahat_curr <- solve(Vinv_curr) %*% crossprod(X_curr, y)
  dcurr <- g2 + ssy - crossprod(ahat_curr, Vinv_curr %*% ahat_curr)
  
  tally <- c(0,0,0)

  for(it in 2:its) {
    
    ## Reversible Jump Step ##
    
    fate <- function(knots = nknot) {
      if((knots == 0) | (knots == 1)) {return(1)} # having 0 or 1 knots must auto defer to birth
      if(knots == max_knot) {sample(2:3, 1)} #at max_knot knot capacity, can only delete or change
      sample(3, 1)
      # 1 = BIRTH
      # 2 = DEATH
      # 3 = CHANGE
    }
    choice <- fate()
    
    ## Eliminating the need for an else statement at the acceptance step ##
    mat_t[i,,] <- mat_t[i-1,,]
    mat_s[i,,] <- mat_s[i-1,,]
    mat_v[i,,] <- mat_v[i-1,,]
    mat_j[i,] <- mat_j[i-1,] #nint
    
    
    samp_j <- function(limit = max_j) {
      sample(1:limit, 1)
    }
    
    samp_vars <- function(deg) {
      p <- ncol(X)
      sample(2:p, deg) %>% sort() #accounting for intercept term, sort the terms
    }
    
    if(choice == 1) { # BIRTH
      # sample j
      # sample j random knots, j random signs, j variables
      # make 1 basis function from this
      # ALMOST the same
      # Prior and proposal will be different
      # evaluate likelihood
      # accept/reject
      j <- samp_j() #nint.cand
      candidate_t <- runif(j, 0, 1) # sample knots for new basis function
      candidate_s <- sample(c(-1,1), j, replace = TRUE)
      vars <- samp_vars(deg = j) 
      Xmat <- X[,vars]
      
      unsign <- t(t(Xmat) - candidate_t)
      basis_mat <- t(candidate_s * t(unsign)) 
      basis_mat <- pos(basis_mat) # matrix of basis functions, same size as X_mat
      # COULD ALSO USE THE spline.basis() FUNCTION
      X_cand <- cbind(X_curr, basis_mat)
      
      Vinv_cand <- crossprod(X_cand) + diag(nknot + 2) / tau2 # +1 for intercept, +1 for presumed birth
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, y)
      dcand <- g2 + ssy - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llik.alpha <- (
        0.5*log(1/tau2) - determinant(Vinv_cand)$mod/2 + determinant(Vinv_curr)$mod/2 + (g1+n/2)*(log(dcurr) - log(dcand))
      ) # calculate the log likelihood ratio
      
      
      
      
      
      
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
        
        V_part <- 0.5 * determinant(Vprime)$mod -  0.5 * determinant(V)$mod
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

  a <- met_gibbs(its = its)
  mat_beta <- a[[1]] #beta values
  mat_sig <- a[[2]] #sig^2 values
  mat_t <- a[[3]]
  mat_s <- a[[4]]
  X_curr <- a[[5]]
  
  list(mat_beta, mat_sig, mat_t, mat_s, X_curr)
}
