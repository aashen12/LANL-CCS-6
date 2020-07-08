
# Full Conditional Sampling With Gibbs
# Andy Shen, Devin Francom
# LANL, CCS-6

### TASK 1 ###
sig_tild <- diag(c(1,4,9))
mu <- mu
S <- diag(3)
M <- c(0,0,0)

## FULL CONDITIONAL DISTRIBUTIONS ##

p_mu_easy <- function(m = M, ST = sig_tild, .S = S, dat = data) {
  sums <- colSums(dat) #update ST
  vals <- rep(NA, 3)
  for(i in seq_len(3)) {
    sig_i <- ST[i,i]
    si <- .S[i,i]
    mi <- m[i]
    mean <- ((sig_i * mi) / (sig_i + n * si)) + 
      ((si * sums[i]) / (sig_i + n * si)) 
    #print(mean)
    var <- (sig_i * si) / (sig_i + n * si)
    vals[i] <- rnorm(1, mean, sqrt(var))
  }
  #return(mean)
  return(vals)
}

p_sig_easy <- function(a = 2, b = 1, dat = data, .mu = mu) {
  centered <- t(t(dat) - .mu) #update .mu
  sums <- colSums(centered^2)
  vals <- rep(NA, 3)
  for(i in seq_len(3)) {
    beta <- 0.5 * sums[i] + b
    vals[i] <- 1/rgamma(1, shape = a + n/2, rate = beta)
  }
  return(vals)
}

## GIBBS SAMPLER ##

gibbs_easy <- function(its = 10000) {
  mat_sig <- matrix(NA, nrow = its, ncol = 3)
  mat_mu <- matrix(NA, nrow = its, ncol = 3)
  mat_sig[1,] <- rep(1, 3)
  mat_mu[1,] <- rep(0, 3)
  for(i in 2:its) {
    mat_sig[i,] <- p_sig_easy(.mu = mat_mu[i-1,])
    mat_mu[i,] <- p_mu_easy(ST = diag(mat_sig[i,]))
  }
  list(mat_mu, mat_sig)
}

#########################################################################################################################

### TASK 2 ###

## FULL CONDITIONALS ##

p_mu <- function(n = nrow(data), ST) {
  covariance <- solve(Sinv + n * solve(ST))
  mean <- covariance %*% (Sinv %*% m + n * solve(ST) %*% colMeans(data))
  rmvnorm(1, mean, covariance)
} #full posterior for mu

p_sig <- function(alpha = 1, beta = 1, ST, sig_sq, .mu) { # Inverse Gamma Prior
  full_cond <- rep(NA, 3)
  
  prior <- (-alpha - 1) * log(sig_sq) + (-beta/sig_sq)
  
  mat <- t(t(data) - .mu) 
  STinv <- tryCatch(solve(ST), error = function(e) browser())
  mat2 <- (mat %*% solve(ST)) * mat
  rsums <- rowSums(mat2)
  
  sum_all <- -0.5 * sum(rsums)
  like <- -0.5 * n * (determinant(ST, logarithm = TRUE)$modulus) + (sum_all)
  
  full_cond <- prior + like
  full_cond
} #full posterior for sigma^2

p_rho <- function(rho, a = 2, b = 3, .mu, ST, n = nrow(data)) { # Beta Prior 
  full_cond <- rep(NA, 3)
  
  rho_trans <- (rho * 0.5) + 0.5
  
  prior <- (a-1) * log(rho_trans) +  (b-1) * log(1 - rho_trans)
  
  mat <- t(t(data) - .mu) 
  mat2 <- (mat %*% solve(ST)) * mat
  rsums <- rowSums(mat2)
  sum_all <- -0.5 * sum(rsums)
  like <- (-n/2) * (determinant(ST, logarithm = TRUE)$modulus) + (sum_all)
  
  full_cond <- prior + like
  full_cond
} #full posterior for rho


## GIBBS SAMPLER ##

cor2cov <- function(cor.mat,vars) {
  mat <- matrix(rep(sqrt(vars), length(vars)), nrow = length(vars))
  cor.mat * mat * t(mat)
}

makeCov <- function(cors,vars) {
  cor.mat <- diag(length(vars))
  cor.mat[upper.tri(cor.mat)] <- cor.mat[lower.tri(cor.mat)] <- cors
  mat <- matrix(rep(sqrt(vars), length(vars)), nrow = length(vars))
  cor.mat * mat * t(mat)
}

met_gibbs <- function(its) {
  sd_sig <- 0.70
  sd_rho <- 0.09
  # tuning sd values
  
  propose_sig <- function(x, mn) {
    log(dtruncnorm(x, a = 0, mean = mn, sd = sd_sig))
  } #truncated normal for proposal
  
  propose_rho <- function(x, mn) {
    log(dtruncnorm(x, a = -1, b = 1, mean = mn, sd = sd_rho))
  }
  
  met_mu <- matrix(NA, its, 3)
  met_mu[1,] <- rep(0, 3)
  
  mat_sig <- matrix(NA, nrow = its, ncol = 3)
  mat_rho <- matrix(NA, nrow = its, ncol = 3)
  mat_sig[1,] <- rep(1, 3)
  mat_rho[1,] <- rep(0, 3) #12, 13, 23
  a_sig <- rep(0, 3)
  a_rho <- rep(0, 3)
  
  for(i in 2:its) {
    #sample mu with mvn
    #sample sigma given mu
    #sample rho given sigma and mu & rtruncnorm between (-1,1)
    
    cov <- makeCov(mat_rho[i-1,], mat_sig[i-1,])
    #browser()
    met_mu[i,] <- p_mu(ST = cov) #Done sampling mu
    
    ### SIGMA ###
    
    for(j in 1:3) {
      candidate_sig <- rtruncnorm(
        1, a = 0, mean = mat_sig[i-1,j], sd = sd_sig
      )
      
      # if(j == 1){candidate_sig = 1}
      # if(j == 2){candidate_sig = 4}
      # if(j == 3){candidate_sig = 9}
      
      cand.sig.vec <- diag(cov)
      cand.sig.vec[j] <- candidate_sig
      cov2 <- makeCov(mat_rho[i-1,], cand.sig.vec)
      
      ratio_sig <- 
        ( p_sig(ST = cov2, sig_sq = candidate_sig, .mu = met_mu[i,]) - 
            propose_sig(candidate_sig, mn = cov[j,j]) ) - 
        ( p_sig(ST = cov, sig_sq = cov[j,j], .mu = met_mu[i,]) - 
            propose_sig(cov[j,j], mn = candidate_sig) )
      
      pmove <- min(0, ratio_sig)
      #if(is.na(pmove)) browser()
      u <- log(runif(1))
      
      if(u < pmove) {
        mat_sig[i,j] <- candidate_sig
        cov2[j,j] <- candidate_sig
        a_sig[j] <- a_sig[j] + 1
        cov <- cov2
      }
      else {
        mat_sig[i,j] <- mat_sig[i-1,j]
      }
    }
    
    ### RHO ###
    
    #browser()
    for(k in 1:3) {
      candidate_rho <- rtruncnorm(
        1, a = -1, b = 1, mean = mat_rho[i-1,k], sd = sd_rho
      )
      
      cor <- cov2cor(cov)
      cand.rho.vec <- cor[upper.tri(cor)]
      cand.rho.vec[k] <- candidate_rho
      
      cov3 <- makeCov(cand.rho.vec, mat_sig[i,])
      
      ratio_cov <- 
        (p_rho(rho = candidate_rho, .mu = met_mu[i,], ST = cov3) -
           propose_rho(candidate_rho, mn = mat_rho[i-1,k])) -
        (p_rho(rho = mat_rho[i-1,k], .mu = met_mu[i,], ST = cov) -
           propose_rho(mat_rho[i-1,k], mn = candidate_rho))
      
      #if(is.na(ratio_cov)) browser()
      
      pmove <- min(0, ratio_cov)
      u <- log(runif(1))
      if(u < pmove) {
        mat_rho[i,k] <- candidate_rho
        a_rho[k] <- a_rho[k] + 1
        cov <- cov3
      }
      else {
        mat_rho[i,k] <- mat_rho[i-1,k]
      }
      #chol(cov)
      #browser()
    }
  }
  list(met_mu, mat_sig, mat_rho, a_sig, a_rho, cov)
}

