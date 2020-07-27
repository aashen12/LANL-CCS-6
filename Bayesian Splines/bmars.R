library(mnormt)
pos<-function(vec){
  (abs(vec)+vec)/2
}
makeBasis<-function(signs,vars,knots,datat,degree=1){
  temp1<-pos(signs*(datat[vars,,drop=F]-knots))^degree # this only works for t(data)...
  if(length(vars)==1){
    return(c(temp1))
  } else{
    temp2<-1
    for(pp in 1:length(vars)){ # faster than apply
      temp2<-temp2*temp1[pp,]
    }
    return(temp2)
  }
}
bass<-function(X,y,max.int=3,max.basis=50,tau2=10^4,nmcmc=10000,g1=0,g2=0,h1=10,h2=10){
  Xt<-t(X)
  n<-length(y)
  p<-ncol(X)
  ssy<-sum(y^2)
  knots<-signs<-vars<-array(dim=c(nmcmc,max.basis,max.int))
  nint<-matrix(nrow=nmcmc,ncol=max.basis) #J
  beta<-matrix(nrow=nmcmc,ncol=max.basis+1) # +1 for intercept
  s2<-lam<-nbasis<-rep(NA,nmcmc)
  nbasis[1]<-0
  s2[1]<-1
  lam[1]<-1
  X.curr<-matrix(rep(1,n))

  Vinv.curr<-crossprod(X.curr)+1/tau2 
  bhat<-solve(Vinv.curr)%*%crossprod(X.curr,y)
  d.curr <- g2 + ssy - crossprod(bhat,Vinv.curr%*%bhat)
  
  count<-c(0,0,0) # count how many times we accept birth, death, change
  beta[1,1]<-bhat
  
  for(i in 2:nmcmc){
    
    ## Reversible jump step
    
    move.type<-sample(c('birth','death','change'),1)
    if(nbasis[i-1]==0)
      move.type<-'birth'
    if(nbasis[i-1]==max.basis)
      move.type<-sample(c('death','change'),1)
    
    # set all of this iterations values to last iteration values...we'll change them if we accpet a move below
    nbasis[i]<-nbasis[i-1]
    nint[i,]<-nint[i-1,]
    knots[i,,]<-knots[i-1,,]
    signs[i,,]<-signs[i-1,,]
    vars[i,,]<-vars[i-1,,]
    
    if(move.type=='birth'){
      nint.cand<-sample(max.int,1) # sample degree of interaction for new basis function
      knots.cand<-runif(nint.cand) # sample knots for new basis function
      signs.cand<-sample(c(-1,1),nint.cand,replace = T) # signs for new basis function
      vars.cand<-sample(p,nint.cand,replace = F) # variables to use in new basis function
      basis.cand<-makeBasis(signs.cand,vars.cand,knots.cand,Xt) # make the new basis function
      X.cand<-cbind(X.curr,basis.cand) # add the new basis function to the basis functions we already have
      
      Vinv.cand<-crossprod(X.cand)+diag(nbasis[i-1]+2)/tau2 # +2: one for intercept and one for birth
      bhat.cand<-solve(Vinv.cand)%*%t(X.cand)%*%y
      d.cand <- g2 + ssy - t(bhat.cand)%*%Vinv.cand%*%bhat.cand
      
      llik.alpha <- .5*log(1/tau2) - determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 + (g1+n/2)*(log(d.curr) - log(d.cand)) # calculate the log likelihood ratio
      lprior.alpha <- ( # log prior ratio
        log(lam[i-1])-log(nbasis[i-1]+1) # nbasis
        + log(1/max.int)  # nint
        + nint.cand*log(1/2)  # signs
        + log(1/choose(p,nint.cand)) # vars
        + log(nbasis[i-1]+1) # ordering
      )
      lprop.alpha <- ( # log proposal ratio
        # probability of going from proposed to current (death)
        log(1/3) + # probability of selection a death step
        + log(1/(nbasis[i-1]+1)) # probability that this basis function is selected to kill
        # probability of going from current to proposed
        - (
            log(1/3) # probability of selecting a birth step
            + log(1/max.int) # probability of nint.cand
            + nint.cand*log(1/2) # probability of the signs we got
            + log(1/choose(p,nint.cand)) # probability of vars.cand
          )
      )
      
      alpha <- llik.alpha + lprior.alpha + lprop.alpha

      if(log(runif(1))<alpha){
        X.curr<-X.cand
        bhat<-bhat.cand
        Vinv.curr<-Vinv.cand
        d.curr<-d.cand
        nbasis[i]<-nbasis[i-1]+1
        nint[i,nbasis[i]]<-nint.cand
        knots[i,nbasis[i],1:nint.cand]<-knots.cand
        signs[i,nbasis[i],1:nint.cand]<-signs.cand
        vars[i,nbasis[i],1:nint.cand]<-vars.cand
        count[1]<-count[1]+1
      }
      
    } else if(move.type=='death'){
      tokill<-sample(nbasis[i-1],1) # which basis function we will delete
      X.cand<-X.curr[,-(tokill+1)] # +1 to skip the intercept
      
      Vinv.cand<-crossprod(X.cand)+diag(nbasis[i-1])/tau2
      bhat.cand<-solve(Vinv.cand)%*%crossprod(X.cand,y)
      d.cand <- g2 + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)
      
      llik.alpha <- -.5*log(1/tau2) - determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 + (g1+n/2)*(log(d.curr) - log(d.cand))
      lprior.alpha <- (
        -log(lam[i-1])+log(nbasis[i-1]) # nbasis
        - log(1/max.int)  # nint
        - nint[i-1,tokill]*log(1/2)  # signs
        - log(1/choose(p,nint[i-1,tokill])) # vars
        - log(nbasis[i-1]) # ordering
      )
      lprop.alpha <- (
        # probability of going from proposed to current (birth)
        log(1/3) # probability of selecting a birth step
        + log(1/max.int) # probability of nint
        + nint[i-1,tokill]*log(1/2) # probability of the signs we got
        + log(1/choose(p,nint[i-1,tokill])) # probability of vars
        # probability of going from current to proposed
        - (
          log(1/3) # probability of selection a death step
          + log(1/nbasis[i-1]) # probability that this basis function is selected to kill
        )
      )
      
      alpha <- llik.alpha + lprior.alpha + lprop.alpha

      if(log(runif(1))<alpha){
        X.curr<-X.cand
        bhat<-bhat.cand
        Vinv.curr<-Vinv.cand
        d.curr<-d.cand
        nbasis[i]<-nbasis[i-1]-1
        nint[i,]<-NA
        knots[i,,]<-signs[i,,]<-vars[i,,]<-NA
        if(nbasis[i]==0){
          nint[i,]<-NA
        } else{
          nint[i,1:nbasis[i]]<-nint[i-1,(1:nbasis[i-1])[-tokill]]
          knots[i,1:nbasis[i],]<-knots[i-1,(1:nbasis[i-1])[-tokill],]
          signs[i,1:nbasis[i],]<-signs[i-1,(1:nbasis[i-1])[-tokill],]
          vars[i,1:nbasis[i],]<-vars[i-1,(1:nbasis[i-1])[-tokill],]
        }
        count[2]<-count[2]+1
      }
      
    } else{
      
      tochange<-sample(nbasis[i-1],1) # which basis function we will change
      tochange2<-sample(nint[i-1,tochange],1) # which element in the basis function tensor product we will change
      knots.cand<-knots[i-1,tochange,1:nint[i-1,tochange]] # copy previous
      knots.cand[tochange2]<-runif(1) # change one element
      signs.cand<-signs[i-1,tochange,1:nint[i-1,tochange]]
      signs.cand[tochange2]<-sample(c(-1,1),1)
      basis<-makeBasis(signs.cand,vars[i-1,tochange,1:nint[i-1,tochange]],knots.cand,Xt)
      X.cand<-X.curr
      X.cand[,tochange+1]<-basis # +1 for intercept
      
      Vinv.cand<-crossprod(X.cand)+diag(nbasis[i-1]+1)/tau2
      bhat.cand<-solve(Vinv.cand)%*%crossprod(X.cand,y)
      d.cand <- g2 + ssy - crossprod(bhat.cand,Vinv.cand%*%bhat.cand)
      
      llik.alpha <- -determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 + (g1+n/2)*(log(d.curr) - log(d.cand))
      
      alpha <- llik.alpha
      
      if(log(runif(1))<alpha){
        X.curr<-X.cand
        bhat<-bhat.cand
        Vinv.curr<-Vinv.cand
        d.curr<-d.cand
        knots[i,tochange,1:nint[i,tochange]]<-knots.cand
        signs[i,tochange,1:nint[i,tochange]]<-signs.cand
        count[3]<-count[3]+1
      }
      
    }
    
    ## Gibbs steps
    
    lam[i]<-rgamma(1,h1+nbasis[i],h2+1)
    S<-solve(crossprod(X.curr)/s2[i-1]+diag(nbasis[i]+1)/tau2)
    beta[i,1:(nbasis[i]+1)]<-rmnorm(1,S%*%t(X.curr)%*%y/s2[i-1],S)
    s2[i]<-1/rgamma(1,n/2+g1,rate=g2+.5*sum((y-X.curr%*%beta[i,1:(nbasis[i]+1)])^2))
  }

  return(list(X=X.curr,b=bhat,count=count,knots=knots,signs=signs,vars=vars,nint=nint,nbasis=nbasis,beta=beta,s2=s2,lam=lam))
}

predict.bass<-function(mod,X,nburn=1000){
  nmcmc<-length(mod$nbasis)
  Xt<-t(X)
  pred<-matrix(nrow=nmcmc-nburn,ncol=nrow(X))
  for(i in (nburn+1):nmcmc){
    B<-matrix(nrow=nrow(X),ncol=mod$nbasis[i]+1)
    B[,1]<-1
    for(j in 1:mod$nbasis[i])
      B[,j+1]<-makeBasis(mod$signs[i,j,1:mod$nint[i,j]],mod$vars[i,j,1:mod$nint[i,j]],mod$knots[i,j,1:mod$nint[i,j]],Xt)
    pred[i-nburn,]<-B%*%mod$beta[i,1:(mod$nbasis[i]+1)]
  }
  return(pred)
}

################################################################################################################
# test it out

set.seed(12)
f <-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}

sigma <- 1 # noise sd
n <- 500 # number of observations
x <- matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y <- rnorm(n,f(x),sigma)

mod<-bass(x,y)
mod$count
plot(mod$nbasis,type='l')
plot(y,mod$X%*%mod$b); abline(a=0,b=1,col=2)
plot(mod$s2[-c(1:1000)],type='l')
plot(mod$lam[-c(1:1000)],type='l')


xtest<-matrix(runif(1000*10),1000,10)
pred<-predict.bass(mod,xtest)
plot(f(xtest),colMeans(pred)); abline(a=0,b=1,col=2)
