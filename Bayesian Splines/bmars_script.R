library(mvtnorm)
library(mnormt)
library(dplyr)

### [a]+ FUNCTION ### 

pos <- function(vec) {
  ((abs(vec) + vec) / 2)
} # set anything less than 0 to 0; creating basis functions [a]+

############################################################################################################

### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(signs, vars, knots, tdat, deg = 1) {
  temp1 <- pos(signs * (tdat[vars,,drop = F] - knots))^deg
  
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

spline.basis2 <- function(signs, vars, knots, tdat, deg = 1) {
  temp1 <- pos(signs * (tdat[vars,,drop = T] - knots))^deg
  
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

############################################################################################################

bmars <- function(X, its, max_knot = 50, max_j = 3, tau2 = 10^4, g1=0, g2=0, h1=10, h2=10) {
  Xt <- t(X)
  n <- length(y)
  p <- ncol(X)
  ssy <- sum(y^2)
  p_t <- function(sig_sq, betahat, X) {
    (-1/(2*sig_sq)) * t(y - (X %*% betahat)) %*% (y - (X %*% betahat))
  } # full conditional for the knot locations, t
  
  p_sig <- function(a = g1, b = g2, X, beta) {
    n <- nrow(X)
    a_term <- a + (n/2)
    b_term <- 0.5 * (2*b + (t(y - (X %*% beta)) %*% (y - (X %*% beta))))
    1 / rgamma(1, shape = a_term, rate = b_term)
  } # full conditional for sigma^2 in the regression equation
  
  p_beta <- function(sig_sq, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau2) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    rmvnorm(1, mean = mu, sigma = sig)
  } # full conditional for the regression coefficients
  
  bhat <- function(sig_sq, X) {
    p = ncol(X) - 1
    sig <- solve( (1/sig_sq) * (t(X) %*% X) + (1/tau2) * diag(p+1) )
    mu <- (1/sig_sq) * sig %*% t(X) %*% y
    mu
  } #marginalized betahat for regression coefficients
  # THIS IS THE GAUSSIAN LIKELIHOOD
    
  ## Initializing Values for RJMCMC ##
  
  mat_t <- array(NA, dim = c(its, max_knot, max_j)) #knots
  mat_s <- array(NA, dim = c(its, max_knot, max_j)) #signs
  mat_v <- array(NA, dim = c(its, max_knot, max_j)) #vars
  
  mat_j <- matrix(NA, its, max_knot) #number of interactions: nint
  mat_beta <- matrix(NA, its, max_knot + 1)
  mat_beta[1,1] <- mean(y)
  mat_sig <- lam <- nknot <- rep(NA, its)
  mat_sig[1] <- 1
  lam[1] <- 1
  
  nknot[1] <- 0 # First iteration must be a birth
  X_curr <- rep(1, n) %>% as.matrix() 
  # first current X-matrix is just an intercept

  Vinv_curr <- crossprod(X_curr) + 1/tau2
  ahat_curr <- solve(Vinv_curr) %*% crossprod(X_curr, y)
  dcurr <- g2 + ssy - crossprod(ahat_curr, Vinv_curr %*% ahat_curr)
  
  count <- c(0,0,0)

  for(i in 2:its) {
    
    ## Reversible Jump Step ##
    
    fate <- function(knots = nknot[i-1]) {
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
    nknot[i] <- nknot[i-1]
    
    
    samp_j <- function(limit = max_j) {
      sample(1:limit, 1)
    }
    
    samp_vars <- function(deg) {
      p <- ncol(X)
      sample(p, deg)
    }
    
    if(choice == 1) { # BIRTH
      j <- sample(max_j, 1) #nint.cand: degree of interaction
      candidate_t <- runif(j, 0, 1) # sample knots for new basis function
      candidate_s <- sample(c(-1,1), j, replace = TRUE)
      var <- samp_vars(deg = j) #vars.cand this is a candidate value
      Xmat <- X[,var]
      
      basis_mat <- spline.basis(signs = candidate_s, vars = var, knots = candidate_t, tdat = Xt)
      X_cand <- cbind(X_curr, basis_mat)
      
      Vinv_cand <- crossprod(X_cand) + diag(nknot[i-1] + 2) / tau2 # +1 for intercept, +1 for presumed birth
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, y)
      dcand <- g2 + ssy - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llik.alpha <- (
        0.5*log(1/tau2) - determinant(Vinv_cand)$mod/2 + determinant(Vinv_curr)$mod/2 + (g1+n/2)*(log(dcurr) - log(dcand))
      ) # calculate the log likelihood ratio
      
      lprior.alpha <- ( # log prior ratio
        log(lam[i-1]) - log(nknot[i-1]+1) # nbasis
        + log(1/max_j)  # nint
        + j * log(1/2)  # signs
        + log(1/choose(p, j)) # vars
        + log(nknot[i-1] + 1) # ordering
      )  
      
      lprop.alpha <- ( # log proposal ratio
        # probability of going from proposed to current (death)
        log(1/3) # probability of selection a death step
        + log(1/(nknot[i-1]+1)) # probability that this basis function is selected to kill
        # probability of going from current to proposed
        - (
          log(1/3) # probability of selecting a birth step
          + log(1/max_j) # probability of j
          + j * log(1/2) # probability of the signs we got
          + log(1/choose(p, j)) # probability of vars.cand
        )
      )
      
      ratio_rj <- llik.alpha + lprior.alpha + lprop.alpha
      #accept_prob <- min(0, ratio_rj)
      
      if(is.na(ratio_rj)) browser()
      if(log(runif(1)) < ratio_rj) { 
        ahat_curr <- ahat_cand
        Vinv_curr <- Vinv_cand
        dcurr <- dcand
        X_curr <- X_cand
        nknot[i] <- nknot[i-1] + 1
        mat_j[i, nknot[i]] <- j
        mat_t[i, nknot[i], 1:j] <- candidate_t
        mat_s[i, nknot[i], 1:j] <- candidate_s
        mat_v[i, nknot[i], 1:j] <- var
        count[1] <- count[1] + 1
      }
    } 
    
    else if(choice == 2) { #DEATH
      # select function to delete
      # prior, proposal different
      pick <- sample(nknot[i-1], 1) 
      X_cand <- X_curr[,-(pick+1)] #ACCOUNTING FOR THE INTERCEPT
      
      Vinv_cand <- crossprod(X_cand) + diag(nknot[i-1]) / tau2
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, y)
      dcand <- g2 + ssy - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llike <- (
        -0.5*log(1/tau2) 
        - determinant(Vinv_cand)$mod + determinant(Vinv_curr)$mod
        + (g1 + n/2) * (log(dcurr) - log(dcand))
      )
      
      lprior <- (
        -log(lam[i-1]) + log(nknot[i-1])
        - log(1/max_j)
        - mat_j[i-1, pick] * log(1/2)
        - log(1/choose(p, mat_j[i-1, pick]))
        - log(nknot[i-1])
      )
      
      lprop <- (
        # prob of going from proposed to current
        log(1/3) #prob of selecting birth
        + log(1/max_j) #prob of j
        + mat_j[i-1, pick] * log(1/2) # probability of signs
        + log(1/choose(p, mat_j[i-1, pick])) #prob of vars
        # prob of going from current to proposed
        - (
          log(1/3)
          + log(1/nknot[i-1])
        )
      )
      
      ratio_rj <- llike + lprior + lprop
      
      #accept_prob <- min(0, ratio_rj)
      
      if(is.na(ratio_rj)) browser()
      
      if(log(runif(1)) < ratio_rj) { 
        X_curr <- X_cand
        ahat_curr <- ahat_cand
        Vinv_curr <- Vinv_cand
        dcurr <- dcand
        nknot[i] <- nknot[i-1] - 1
        mat_j[i,] <- NA
        mat_t[i,,] <- mat_s[i,,] <- mat_v[i,,] <- NA
        if(nknot[i] == 0) {
          mat_j[i,] <- NA
        } else {
          mat_j[i,(1:nknot[i])] <- mat_j[i-1,(1:nknot[i-1])[-pick]]
          mat_t[i,(1:nknot[i]),] <- mat_t[i-1,(1:(nknot[i-1]))[-pick],]
          mat_s[i,(1:nknot[i]),] <- mat_s[i-1,(1:(nknot[i-1]))[-pick],]
          mat_v[i,(1:nknot[i]),] <- mat_v[i-1,(1:(nknot[i-1]))[-pick],]
        }
        count[2] <- count[2] + 1
      } 
    } 
    
    else { #CHANGE
     
      tochange <- sample(nknot[i-1], 1)
      tochange2 <- sample(mat_j[i-1,tochange], 1)
      candidate_t <- mat_t[i-1, tochange, 1:mat_j[i-1, tochange]] #copy previous
      candidate_t[tochange2] <- runif(1) #change 1 element
      candidate_s <- mat_s[i-1, tochange, 1:mat_j[i-1, tochange]]
      candidate_s[tochange2] <- sample(c(-1, 1), 1) #change 1 element
      # same idea
      
      basis <- spline.basis(candidate_s, mat_v[i-1,tochange,1:mat_j[i-1,tochange]], candidate_t, Xt)
      
      X_cand <- X_curr
      X_cand[,tochange+1] <- basis #+1 for int
      
      Vinv_cand <- crossprod(X_cand) + diag(nknot[i-1]+1)/tau2
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, y)
      dcand <- g2 + ssy - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llik <- (
        -determinant(Vinv_cand)$mod/2 
        + determinant(Vinv_curr)$mod/2 
        + (g1+n/2) * (log(dcurr) - log(dcand))
      )
      
      #acc <- min(0, llik)
      
      if(log(runif(1)) < llik) {
        X_curr <- X_cand
        ahat_curr <- ahat_cand
        Vinv_curr <- Vinv_cand
        dcurr <- dcand
        mat_t[i,tochange,1:mat_j[i,tochange]] <- candidate_t
        mat_s[i,tochange,1:mat_j[i,tochange]] <- candidate_s
        count[3] <- count[3] + 1
      }
    }
    lam[i] <- rgamma(1, h1 + nknot[i], h2 + 1)
    S <- solve(crossprod(X_curr) / mat_sig[i-1] + diag(nknot[i] + 1)/tau2)
    mat_beta[i,1:(nknot[i]+1)] <- rmnorm(1, S %*% t(X_curr) %*% y/mat_sig[i-1], S)
    mat_sig[i] <- 1/rgamma(
      1, 
      shape = (n/2) + g1,
      rate = g2 + 0.5 * sum((y-X_curr %*% mat_beta[i,1:(nknot[i]+1)])^2)
    )
  }
  
  names(count) <- c("birth", "death", "change")
  
  list(X = X_curr, beta = ahat_curr, 
       count = count, knots = mat_t, signs = mat_s, 
       vars = mat_v, int = mat_j, nknot = nknot, mat_beta = mat_beta, 
       mat_sig = mat_sig, lam = lam)

}

############################################################################################################

### PREDICTION FUNCTION ###

predict.bmars <- function(mod, X, nburn = 1000) {
  nmcmc <- length(mod$nknot)
  Xt <- t(X)
  pred <- matrix(NA, nrow = nmcmc-nburn, ncol = nrow(X))
  
  for(i in (nburn+1):nmcmc) {
    B <- matrix(NA, nrow(X), mod$nknot[i] + 1)
    B[,1] <- 1
    for(j in (1:mod$nknot[i])) {
      #browser()
      B[,j+1] <- spline.basis(signs = mod$signs[i,j,1:mod$int[i,j]],
                              vars = mod$vars[i,j,1:mod$int[i,j]],
                              knots = mod$knots[i,j,1:mod$int[i,j]],
                              tdat = Xt)
      
    }
    #browser()
    pred[i-nburn,] <- B %*% mod$mat_beta[i,1:(mod$nknot[i]+1)]
  }
  
  pred
}
