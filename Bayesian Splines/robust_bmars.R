# Function to find the optimal adaptive regression spline for multivariate data assuming non-Gaussian likelihoods
# Andy Shen, Devin Francom
# Statistical Sciences Group (CCS-6): Los Alamos National Laboratory

############################################################################################################

library(mvtnorm)
library(mnormt)
library(dplyr)
library(LaplacesDemon)

### [a]+ FUNCTION ### 

pos <- function(vec) {
  ((abs(vec) + vec) / 2)
} # set anything less than 0 to 0; creating basis functions [a]+

############################################################################################################

### SPLINE BASIS FUNCTION ###

# Generate a matrix of basis functions given number of knots and vector of knots

spline.basis <- function(signs, vars, knots, tdat, deg = 1) {
  temp1 <- pos(signs * (tdat[vars,,drop = F] - knots))^deg #transformation
  
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
############################################################################################################

### BMARS ALGORITHM ###

bmars <- function(X, its, max_knot=50, max_j=3, tau2=10^4, g1=0, g2=0, h1=10, h2=10, nu=10) {
  Xt <- t(X)
  n <- length(y)
  p <- ncol(X)
  ssy <- sum(y^2)
  
  ## Initializing Values for RJMCMC ##
  
  mat_t <- array(NA, dim = c(its, max_knot, max_j)) #knots
  mat_s <- array(NA, dim = c(its, max_knot, max_j)) #signs
  mat_v <- array(NA, dim = c(its, max_knot, max_j)) #vars
  
  mat_w <- array(NA, dim = c(n, n, its)) #an nxn matrix replicated for each iteration
  
  mat_j <- matrix(NA, its, max_knot) #number of interactions: nint
  mat_beta <- matrix(NA, its, max_knot + 1) #+1 for intercept
  mat_beta[1,1] <- mean(y)
  mat_sig <- lam <- nknot <- rep(NA, its) #sigma^2, lambda, number of knots
  mat_sig[1] <- 1
  lam[1] <- 1
  
  mat_w[,,1] <- diag(rep(1, n))
  
  nknot[1] <- 0 # First iteration must be a birth
  X_curr <- rep(1, n) %>% as.matrix() # first current X-matrix is just an intercept
  
  Wcurr <- mat_w[,,1]
  
  ### Likelihood in Denison, et. al ###
  Vinv_curr <- crossprod(X_curr, Wcurr %*% X_curr) + 1/tau2 
  ahat_curr <- solve(Vinv_curr) %*% crossprod(X_curr, Wcurr %*% y)
  dcurr <- g2 + ssy - crossprod(ahat_curr, Vinv_curr %*% ahat_curr)
  
  count <- c(0,0,0) #track number of birth, death and change moves
  
  for(i in 2:its) {
    
    ## Reversible Jump Step ##
    
    # Uniformly decide whether to add, delete or change a basis function
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
    mat_t[i,,] <- mat_t[i-1,,] #knot locations
    mat_s[i,,] <- mat_s[i-1,,] #signs
    mat_v[i,,] <- mat_v[i-1,,] #vars
    mat_j[i,] <- mat_j[i-1,] #nint
    nknot[i] <- nknot[i-1]
    
    mat_w[,,i] <- mat_w[,,i-1]
    
    # samp_j <- function(limit = max_j) {
    #   sample(1:limit, 1)
    # }
    
    samp_vars <- function(deg) {
      p <- ncol(X)
      sample(p, deg)
    } #based on j, sample the actual columns of your X matrix that you will consider
    
    if(choice == 1) { # BIRTH
      j <- sample(max_j, 1) #nint.cand: degree of interaction for new basis function
      candidate_t <- runif(j, 0, 1) # sample knot locations for new basis function
      candidate_s <- sample(c(-1,1), j, replace = TRUE) #signs for new basis functions
      
      candidate_w <- 1/rchisq(n, nu+1) #sample candidate w from inv Chi Sq
      # candidate_w <- rinvchisq(
      #   n,
      #   df = nu + 1,
      #   scale = ((nu * mat_sig[i-1]) + (y - X %*% mat_beta[i-1,])^2) / (nu + 1)
      # )
      Wcand <- diag(1/candidate_w) #matrix is 1/candidate
      
      var <- samp_vars(deg = j) #vars.cand: this is a candidate value for the columns of X to extract. There are j of them
      basis_mat <- spline.basis(signs = candidate_s, vars = var, knots = candidate_t, tdat = Xt) #candidate basis function
      X_cand <- cbind(X_curr, basis_mat)
      
      ### BASED ON DENISON, ET. AL ###
      
      Vinv_cand <- crossprod(X_cand, Wcand %*% X_cand) + diag(nknot[i-1] + 2) / tau2 # +2: 1 for intercept, 1 for birth
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, Wcand %*% y)
      dcand <- g2 + crossprod(y, Wcand%*%y) - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      #update the likelihood based on t-errors
      
      llik.alpha <- (
        0.5*log(1/tau2) # simplifying the tau2*I fraction
        - determinant(Vinv_cand)$mod/2 + determinant(Vinv_curr)$mod/2 
        + (g1+n/2)*(log(dcurr) - log(dcand))
      ) # calculate the log likelihood ratio
      
      lprior.alpha <- ( # log prior ratio
        log(lam[i-1]) - log(nknot[i-1]+1) # nbasis: lambda/(M+1)
        + log(1/max_j)  # nint
        + j * log(1/2)  # signs
        + log(1/choose(p, j)) # vars
        + log(nknot[i-1] + 1) # ordering
      )  
      
      lprop.alpha <- ( # log proposal ratio
        # probability of going from proposed to current (death)
        log(1/3) # probability of selecting a death step
        + log(1/(nknot[i-1]+1)) # probability that this basis function is selected to delete
        # probability of going from current to proposed
        - (
          log(1/3) # probability of selecting a birth step
          + log(1/max_j) # probability of j
          + log(1/choose(p, j)) # probability of vars.cand
          + j * log(1/2) # probability of the signs we got
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
        Wcurr <- Wcand #accept the W matrix
        nknot[i] <- nknot[i-1] + 1 #tracking the knots at every iteration
        mat_j[i, nknot[i]] <- j #the candidate degree of interaction is accepted
        mat_t[i, nknot[i], 1:j] <- candidate_t
        mat_s[i, nknot[i], 1:j] <- candidate_s
        mat_v[i, nknot[i], 1:j] <- var
        mat_w[,,i] <- Wcand #update W matrix
        count[1] <- count[1] + 1
      }
    } 
    
    else if(choice == 2) { #DEATH
      pick <- sample(nknot[i-1], 1) # select a knot to delete from the current number
      X_cand <- X_curr[,-(pick+1)] # ACCOUNTING FOR THE INTERCEPT as pick is based on number of knots. 
      #Wcand <- Wcurr[,-pick]
      # ncol(X) = nknot[i-1] + 1
      # HOW DO YOU DELETE A BASIS FUNCTION FROM W
      
      ### BASED ON DENISON, ET. AL ###
      
      Vinv_cand <- crossprod(X_cand, Wcurr %*% X_cand) + diag(nknot[i-1]+1-1) / tau2 #+1 int, -1 death
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, Wcurr %*% y)
      dcand <- g2 + crossprod(y, Wcurr%*%y) - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llike <- (
        -0.5*log(1/tau2) 
        - determinant(Vinv_cand)$mod/2 + determinant(Vinv_curr)$mod/2
        + (g1 + n/2) * (log(dcurr) - log(dcand))
      )
      
      lprior <- (
        -log(lam[i-1]) + log(nknot[i-1])
        - log(1/max_j)
        - mat_j[i-1, pick] * log(1/2) # J_m => previous row (iteration), column of desired death as ncol = max_knot
        - log(1/choose(p, mat_j[i-1, pick]))
        - log(nknot[i-1]) # M
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
        # HOW IS THE W MATRIX UPDATED IN THE DEATH/CHANGE STEPS?
        X_curr <- X_cand
        ahat_curr <- ahat_cand
        Vinv_curr <- Vinv_cand
        dcurr <- dcand
        nknot[i] <- nknot[i-1] - 1
        mat_j[i,] <- NA #no degree of interaction for a death
        mat_t[i,,] <- mat_s[i,,] <- mat_v[i,,] <- NA #no signs, knot locations, or variables involved
        if(nknot[i] == 0) {
          mat_j[i,] <- NA # if you now have no knots, degree of interaction is NA
        } else {
          mat_j[i,(1:nknot[i])] <- mat_j[i-1,(1:nknot[i-1])[-pick]] #new degree of interaction is the old without the selected one
          mat_t[i,(1:nknot[i]),] <- mat_t[i-1,(1:(nknot[i-1]))[-pick],] #etc
          mat_s[i,(1:nknot[i]),] <- mat_s[i-1,(1:(nknot[i-1]))[-pick],] #etc
          mat_v[i,(1:nknot[i]),] <- mat_v[i-1,(1:(nknot[i-1]))[-pick],] #etc
          # this is updated for every sheet in each array since degree of interaction is irrelevant
        }
        count[2] <- count[2] + 1
      } 
    } 
    
    else { #CHANGE
      # HOW DO YOU CHANGE A BASIS FUNCTION FROM W?
      tochange <- sample(nknot[i-1], 1) #sample a knot number to change 
      tochange2 <- sample(mat_j[i-1,tochange], 1) #sample a number from the...
      #degree of interaction corresponding to this knot number,
      # effectively a candidate j. Sampling a degree of interaction to change
      candidate_t <- mat_t[i-1, tochange, 1:mat_j[i-1, tochange]] #copy previous
      candidate_t[tochange2] <- runif(1) #change 1 element, tochange2 %in% (1:mat_j[i-1, tochange])
      candidate_s <- mat_s[i-1, tochange, 1:mat_j[i-1, tochange]] #same idea
      candidate_s[tochange2] <- sample(c(-1, 1), 1) #change 1 element
      # no candidate_vars?
      basis <- spline.basis(candidate_s, mat_v[i-1,tochange,1:mat_j[i-1,tochange]], candidate_t, Xt)
      
      X_cand <- X_curr
      X_cand[,tochange+1] <- basis # +1 for int. updating the knot column with the newly proposed basis function
      
      ### BASED ON DENISON, ET. AL ###
      
      Vinv_cand <- crossprod(X_cand, Wcurr %*% Xcand) + diag(nknot[i-1]+1)/tau2
      ahat_cand <- solve(Vinv_cand) %*% crossprod(X_cand, Wcurr %*% y)
      dcand <- g2 + crossprod(y, Wcurr%*%y) - crossprod(ahat_cand, Vinv_cand %*% ahat_cand)
      
      llik <- (
        -determinant(Vinv_cand)$mod/2 
        + determinant(Vinv_curr)$mod/2 
        + (g1+n/2) * (log(dcurr) - log(dcand))
      )
      
      #acc <- min(0, llik)
      if(is.na(llik)) browser()
      if(log(runif(1)) < llik) {
        #nknot[i] = nknot[i-1] has already been done way up above
        X_curr <- X_cand
        ahat_curr <- ahat_cand
        Vinv_curr <- Vinv_cand
        dcurr <- dcand
        mat_t[i,tochange,1:mat_j[i,tochange]] <- candidate_t
        mat_s[i,tochange,1:mat_j[i,tochange]] <- candidate_s
        count[3] <- count[3] + 1
      }
    }
    
    ### GIBBS SAMPLING STEPS FOR CONSTANTS ###  
    
    lam[i] <- rgamma(1, h1 + nknot[i], h2 + 1)
    
    S <- solve(crossprod(X_curr) / mat_sig[i-1] + diag(nknot[i] + 1)/tau2)
    mat_beta[i,1:(nknot[i]+1)] <- (
      solve(t(X_curr) %*% Wcurr %*% X_curr) %*% t(X_curr) %*% Wcurr %*% y
    )
    
    mat_sig[i] <- (1/n)*crossprod(
      (y- X_curr %*% mat_beta[i,1:(nknot[i]+1)]),
      Wcurr %*% (y- X_curr %*% mat_beta[i,1:(nknot[i]+1)])
    )
  }
  
  names(count) <- c("birth", "death", "change")
  
  list(X = X_curr, beta = ahat_curr, 
       count = count, knots = mat_t, signs = mat_s, 
       vars = mat_v, int = mat_j, nknot = nknot, mat_beta = mat_beta, 
       mat_sig = mat_sig, lam = lam, W = Wcurr, mat_w = mat_w)
  
}
############################################################################################################
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

############################################################################################################
############################################################################################################