# Function to find the optimal adaptive regression spline for multivariate data assuming Student's t-distributed likelihoods
# Andy Shen, Devin Francom
# Statistical Sciences Group (CCS-6): Los Alamos National Laboratory

############################################################################################################

library(mvtnorm)
library(mnormt)
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
  #vectorized, hockey sticks
  
  if(length(vars) == 1) {
    return(c(temp1)) #if everything scalar, return vector
  } else {
    temp2 <- 1
    for(pp in 1:length(vars)) { #if knot, sign, vars are vectors
      temp2 <- temp2 * temp1[pp,] #multivariate
    }
    return(temp2) #apply(temp1, 1, prod)
  }
  
}

############################################################################################################
############################################################################################################

### BMARS ALGORITHM ###

bmars <- function(X, its, max_knot=50, max_j=3, tau2=10^4, g1=0, g2=0, h1=10, h2=10, nu=10, verbose = FALSE) {
  Xt <- t(X)
  n <- length(y)
  p <- ncol(X)
  ssy <- sum(y^2)
  
  ## Initializing Values for RJMCMC ##
  
  mat_t <- array(NA, dim = c(its, max_knot, max_j)) #knots
  mat_s <- array(NA, dim = c(its, max_knot, max_j)) #signs
  mat_v <- array(NA, dim = c(its, max_knot, max_j)) #vars
  
  mat_u <- mat_w <- matrix(NA, its, n) #holds diagonals of the W matrix
  
  mat_j <- matrix(NA, its, max_knot) #number of interactions: nint
  mat_beta <- matrix(NA, its, max_knot + 1) #+1 for intercept
  mat_beta[1,1] <- runif(1) #mean(y)
  mat_sig <- mat_tau2 <- lam <- nknot <- a2 <- rep(NA, its) #sigma^2, lambda, number of knots
  mat_tau2[1] <- runif(1) #1
  lam[1] <- runif(1,0,5) #1
  a2[1] <- log(runif(1))^2
  mat_u[1,] <- rep(runif(1),n) #rep(1, n)
  
  mat_w[1,] <- a2[1] * mat_u[1,]
  Wcurr <- diag(mat_w[1,])
  
  mat_sig[1] <- mat_tau2[1] * a2[1]
  nknot[1] <- 0 # First iteration must be a birth
  X_curr <- as.matrix(rep(1,n)) # first current X-matrix is just an intercept
  
  
  # Wcurr inverse is diag(1/diag(Wcurr))
  
  ### Likelihood in Denison, et. al ###
  
  #updated for robust t
  # Vinv_curr <- crossprod(X_curr, Wcurr %*% X_curr) + 1/tau2
  # ahat_curr <- solve(Vinv_curr) %*% crossprod(X_curr, Wcurr %*% y)
  # dcurr <- g2 + crossprod(y, Wcurr%*%y) - crossprod(ahat_curr, Vinv_curr %*% ahat_curr)
  Hinv_curr <- solve(crossprod(X_curr, diag(1/diag(Wcurr))%*%X_curr) + 1/tau2 * diag(ncol(X_curr)))
  bhat_curr <- Hinv_curr %*% t(X_curr) %*% diag(1/diag(Wcurr)) %*% y
  #browser()
  #dcurr <- g2 + ssy - crossprod(bhat_curr, solve(Hinv_curr)%*%bhat_curr)
  #ssy
  count <- c(0,0,0) #track number of birth, death and change moves
  
  for(i in 2:its) {
    
    ## Reversible Jump Step ##
    # Uniformly decide whether to add, delete or change a basis function
    #if(nknot[i-1] == max_knot) browser()
    fate <- function(knots = nknot[i-1]) {
      if((knots == 0)) {return(1)} # having 0 or 1 knots must auto defer to birth
      if(knots == max_knot) {return(sample(2:3, 1))} #at max_knot knot capacity, can only delete or change
      else{return(sample(3, 1))}
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
    
    samp_vars <- function(deg) {
      p <- ncol(X)
      sample(p, deg)
    } #based on j, sample the actual columns of your X matrix that you will consider
    
    if(choice == 1) { # BIRTH
      j <- sample(max_j, 1) #nint.cand: degree of interaction for new basis function
      candidate_t <- runif(j, 0, 1) # sample knot locations for new basis function
      candidate_s <- sample(c(-1,1), j, replace = TRUE) #signs for new basis functions
      var <- samp_vars(deg = j) #vars.cand: this is a candidate value for the columns of X to extract. There are j of them
      basis_mat <- spline.basis(signs = candidate_s, vars = var, knots = candidate_t, tdat = Xt) #candidate basis function
      X_cand <- cbind(X_curr, basis_mat)
      
      Hinv_cand <- solve(
        crossprod(X_cand, diag(1/diag(Wcurr))%*%X_cand) + 1/tau2 * diag(nknot[i-1]+1+1)
      )
      bhat_cand <- Hinv_cand %*% t(X_cand) %*% diag(1/diag(Wcurr)) %*% y
      #dcand <- g2 + ssy - crossprod(bhat_cand, solve(Hinv_cand)%*%bhat_cand)
      
      llik.alpha <- (
        0.5*log(1/tau2) # simplifying the tau2*I fraction
        + determinant(Hinv_cand)$mod - determinant(Hinv_curr)$mod
        + 0.5*(
          crossprod(bhat_cand, solve(Hinv_cand)%*% bhat_cand) 
          - crossprod(bhat_curr, solve(Hinv_curr)%*% bhat_curr)
        )
        #+ (g1+n/2)*(log(dcurr) - log(dcand))
      ) # calculate the log likelihood ratio
      
      #if(is.na(llik.alpha)) browser()
      
      # llik.alpha <- 
      # (-((nknot[i-1]+1)/2) * log(tau2) - 0.5*(crossprod(y, diag(1/diag(Wcurr))%*%y) - crossprod(bhat_cand,solve(Hinv_cand)%*%bhat_cand)))
      # - (-((nknot[i-1])/2) * log(tau2) - 0.5*(crossprod(y, diag(1/diag(Wcurr))%*%y) - crossprod(bhat,solve(Hinv)%*%bhat)))
      
      lprior.alpha <- ( # log prior ratio
        log(lam[i-1]) - log(nknot[i-1]+1) # nbasis: lambda/(M+1)
        + j * log(1/2)  # signs
        + log(1/choose(p, j)) # vars
        + log(1/max_j)  # nint
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
      if(is.na(ratio_rj)) browser()
      if(log(runif(1)) < ratio_rj) { 
        Hinv_curr <- Hinv_cand
        bhat_curr <- bhat_cand
        #dcurr <- dcand
        X_curr <- X_cand
        nknot[i] <- nknot[i-1] + 1 #tracking the knots at every iteration
        mat_j[i, nknot[i]] <- j #the candidate degree of interaction is accepted
        mat_t[i, nknot[i], 1:j] <- candidate_t
        mat_s[i, nknot[i], 1:j] <- candidate_s
        mat_v[i, nknot[i], 1:j] <- var
        count[1] <- count[1] + 1
      }
    } 
    
    else if(choice == 2) { #DEATH
      pick <- sample(nknot[i-1], 1) # select a knot to delete from the current number
      X_cand <- X_curr[,-(pick+1)] # ACCOUNTING FOR THE INTERCEPT as pick is based on number of knots. 
      
      ### BASED ON DENISON, ET. AL ###
      
      Hinv_cand <- solve(
        crossprod(X_cand, diag(1/diag(Wcurr))%*%X_cand) + 1/tau2* diag(nknot[i-1]-1+1)
      )# * diag(ncol(X_cand)))
      bhat_cand <- Hinv_cand %*% t(X_cand) %*% diag(1/diag(Wcurr)) %*% y
      
      llike <- (
        - 0.5*log(1/tau2) # simplifying the tau2*I fraction
        + 0.5*determinant(Hinv_cand)$mod - 0.5*determinant(Hinv_curr)$mod
        + 0.5*(
          crossprod(bhat_cand, solve(Hinv_cand)%*% bhat_cand) 
          - crossprod(bhat_curr, solve(Hinv_curr)%*% bhat_curr)
        )
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
        X_curr <- X_cand
        bhat_curr <- bhat_cand
        Hinv_curr <- Hinv_cand
        #dcurr <- dcand
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
      
      Hinv_cand <- solve(
        crossprod(X_cand, diag(1/diag(Wcurr))%*%X_cand) + 1/tau2* diag(nknot[i-1]+1)
      )
      bhat_cand <- Hinv_cand %*% t(X_cand) %*% diag(1/diag(Wcurr)) %*% y
      
      llik <- (
        + 0.5*determinant(Hinv_cand)$mod - 0.5*determinant(Hinv_curr)$mod
        + 0.5*(
          crossprod(bhat_cand, solve(Hinv_cand)%*% bhat_cand) 
          - crossprod(bhat_curr, solve(Hinv_curr)%*% bhat_curr)
        )
      )
      
      if(is.na(llik)) browser()
      # ACCEPTANCE
      if(log(runif(1)) < llik) {
        #nknot[i] = nknot[i-1] has already been done way up above
        X_curr <- X_cand
        Hinv_curr<-Hinv_cand
        bhat_curr<-bhat_cand
        mat_t[i,tochange,1:mat_j[i,tochange]] <- candidate_t
        mat_s[i,tochange,1:mat_j[i,tochange]] <- candidate_s
        count[3] <- count[3] + 1
      }
    }
    
    ### GIBBS SAMPLING STEPS FOR CONSTANTS ###  
    
    lam[i] <- rgamma(1, h1 + nknot[i], h2 + 1) #pull number of basis down

    
    mat_beta[i,1:(nknot[i]+1)] <- rmnorm(
      1,
      Hinv_curr %*% crossprod(X_curr, diag(1/diag(Wcurr))%*%y),
      Hinv_curr
    )
    #browser()
    mat_u[i,] <- rinvchisq(
      n,
      nu+1,
      ((nu*mat_tau2[i-1]) + ( (y-(as.matrix(X_curr) %*% as.matrix(mat_beta[i,(1:(nknot[i]+1))]) ))/a2[i-1] )^2) / (nu+1)
    )
    
    # mat_u[i,] <- 1/rgamma(
    #   n,
    #   shape = (nu+1)/2,
    #   rate = ( nu*mat_tau2[i-1] + (y - (X_curr%*%bhat_curr))^2 )/2
    # ) #full conditional of V_i
    # #mat_u[i,]<-1
    
    mat_tau2[i] <- rgamma(
      1,
      (n*nu/2),
      ((nu/2)*sum(1/mat_u[i,]))
    )
    
    a2[i] <- rinvchisq(
      1,
      n,
      (1/n)*sum( (y - (as.matrix(X_curr) %*% as.matrix(mat_beta[i,(1:(nknot[i]+1))]) ))^2 / mat_u[i,] )
    )
    
    if(any(mat_u[i,] > 1^100)) {
      mat_u[i,] <- mat_u[i,]/(10^100)
      a2[i] <- a2[i] * 10^100
      mat_tau2[i] <- mat_tau2[i] / 10^100
    }
    
    if(a2[i] < 10^-18) { #scaling
      a2[i] <- a2[i] * 10^6
      mat_u[i,] <- mat_u[i,] / 10^6
      mat_tau2[i] <- mat_tau2[i] / 10^6
    } else if(a2[i] > 10^18) { #scaling
      a2[i] <- a2[i] / 10^6
      mat_u[i,] <- mat_u[i,] * 10^6
      mat_tau2[i] <- mat_tau2[i] * 10^6
    } 
  
    mat_sig[i] <- a2[i] * mat_tau2[i]
    mat_w[i,] <- a2[i] * mat_u[i,]  #v_i
    
    Wcurr <- diag(mat_w[i,])
    
    Hinv_curr <- solve(
      crossprod(X_curr, diag(1/diag(Wcurr))%*%X_curr) + 1/tau2 * diag(nknot[i]+1)
    ) #Hinv_curr must be updated if Wcurr is updated!!!
    
    if(verbose == TRUE) {
      if(i %% 1000 == 0) {
        cat("Iteration number", i, "\n")
      }
    }
    i <- i + 1
  }
  
  names(count) <- c("birth", "death", "change")
  
  list(X = X_curr, beta = bhat_curr, 
       count = count, knots = mat_t, signs = mat_s, 
       vars = mat_v, int = mat_j, nknot = nknot, mat_beta = mat_beta, 
       mat_tau2 = mat_tau2, lam = lam, W = Wcurr, mat_u = mat_u, mat_w = mat_w,
       mat_sig = mat_sig, alpha=a2)
}
############################################################################################################
############################################################################################################
