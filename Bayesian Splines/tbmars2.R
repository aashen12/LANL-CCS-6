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
tbass<-function(X,y,max.int=3,max.basis=50,tau2=10^4,nu=30,nmcmc=10000,g1=0,g2=0,h1=10,h2=10){
  Xt<-t(X)
  n<-length(y)
  p<-ncol(X)
  ssy<-sum(y^2)
  knots<-signs<-vars<-array(dim=c(nmcmc,max.basis,max.int))
  nint<-matrix(nrow=nmcmc,ncol=max.basis) #J
  beta<-matrix(nrow=nmcmc,ncol=max.basis+1) # +1 for intercept
  s2<-lam<-nbasis<-rep(NA,nmcmc)
  v<-matrix(nrow=nmcmc,ncol=n)
  v[1,]<-1
  nbasis[1]<-0
  s2[1]<-1
  lam[1]<-1
  X.curr<-matrix(rep(1,n))
  
  Vinv.curr<-t(X.curr)%*%diag(v[1,])%*%X.curr/s2[1]+1/tau2 
  bhat<-solve(Vinv.curr)%*%t(X.curr)%*%diag(v[1,])%*%y/s2[1]
  d.curr <- t(y)%*%diag(v[1,])%*%y/s2[1] - crossprod(bhat,Vinv.curr%*%bhat)
  
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
      
      Vinv.cand<-t(X.cand)%*%diag(v[i-1,])%*%X.cand/s2[i-1]+diag(nbasis[i-1]+2)/tau2 # +2: one for intercept and one for birth
      bhat.cand<-solve(Vinv.cand)%*%t(X.cand)%*%diag(v[i-1,])%*%y/s2[i-1]
      d.cand <- - t(bhat.cand)%*%Vinv.cand%*%bhat.cand
      
      llik.alpha <- .5*log(1/tau2) - determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 - .5*d.cand + .5*d.curr # calculate the log likelihood ratio
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
      #browser()
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
      
      Vinv.cand<-t(X.cand)%*%diag(v[i-1,])%*%X.cand/s2[i-1]+diag(nbasis[i-1])/tau2
      bhat.cand<-solve(Vinv.cand)%*%t(X.cand)%*%diag(v[i-1,])%*%y/s2[i-1]
      d.cand <- - t(bhat.cand)%*%Vinv.cand%*%bhat.cand
      
      llik.alpha <- -.5*log(1/tau2) - determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 -.5*d.cand + .5*d.curr
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
      
      Vinv.cand<-t(X.cand)%*%diag(v[i-1,])%*%X.cand/s2[i-1]+diag(nbasis[i-1]+1)/tau2
      bhat.cand<-solve(Vinv.cand)%*%t(X.cand)%*%diag(v[i-1,])%*%y/s2[i-1]
      d.cand <- - t(bhat.cand)%*%Vinv.cand%*%bhat.cand
      
      llik.alpha <- -determinant(Vinv.cand)$mod/2 + determinant(Vinv.curr)$mod/2 - .5*d.cand + .5*d.curr
      
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
    S<-solve(t(X.curr)%*%diag(v[i-1,])%*%X.curr/s2[i-1]+diag(nbasis[i]+1)/tau2)
    beta[i,1:(nbasis[i]+1)]<-rmnorm(1,S%*%t(X.curr)%*%diag(v[i-1,])%*%y/s2[i-1],S)
    v[i,]<-rgamma(n,(nu+1)/2,nu/2 + .5/s2[i-1]*(y-X.curr%*%t(beta[i,1:(nbasis[i]+1),drop=F]))^2)
    s2[i]<-1/rgamma(1,n/2+g1,rate=.5*sum(v[i,]*(y-X.curr%*%t(beta[i,1:(nbasis[i]+1),drop=F]))^2))
    
    Vinv.curr<-t(X.curr)%*%diag(v[i,])%*%X.curr/s2[i]+diag(nbasis[i]+1)/tau2
    bhat.curr<-solve(Vinv.curr)%*%t(X.curr)%*%diag(v[i,])%*%y/s2[i]
    d.curr <- - t(bhat.curr)%*%Vinv.curr%*%bhat.curr
    
    if(i%%1000==0)
      cat(timestamp(quiet = T),' nbasis: ',nbasis[i],'\n')
    
  }
  
  return(list(X=X.curr,b=bhat,count=count,knots=knots,signs=signs,vars=vars,nint=nint,nbasis=nbasis,beta=beta,s2=s2,lam=lam,v=v))
}

predict.tbass<-function(mod,X,nburn=1000){
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

f <-function(x){
  10*sin(pi*x[,1]*x[,2])+20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
}
set.seed(1)
sigma <- 1 # noise sd
n <- 500 # number of observations
x <- matrix(runif(n*10),n,10) #10 variables, only first 5 matter
y <- rnorm(n,f(x),sigma)
ind<-sample(n,size=10)
y[ind]<-rnorm(10,f(x[ind,]),10)
col<-rep(1,n)
col[ind]<-2

mod<-tbass(x,y,nu=10)
mod$count
plot(mod$nbasis,type='l')
plot(y,mod$X%*%mod$b,col=col); abline(a=0,b=1,col=2)
plot(sqrt(mod$s2[-c(1:1000)]),type='l')
matplot(sqrt(1/mod$v[seq(1000,10000,100),]),type='l')
plot(mod$lam[-c(1:1000)],type='l')

hist(colMeans(1/mod$v[-c(1:5000),]),breaks=50)
plot(colMeans(1/mod$v[-c(1:5000),]),(y-mod$X%*%mod$b)^2,col=col)

qq<-quantile(colMeans(sqrt(1/mod$v[-c(1:5000),])),probs = .99)
plot(y,mod$X%*%mod$b,col=col); abline(a=0,b=1,col=2)
pred.out<-which(colMeans(sqrt(1/mod$v[-c(1:5000),])) > qq)
points(y[pred.out],(mod$X%*%mod$b)[pred.out],col=3,pch='+')


xtest<-matrix(runif(1000*10),1000,10)
pred<-predict.tbass(mod,xtest)
plot(f(xtest),colMeans(pred)); abline(a=0,b=1,col=2)
var(f(xtest)-colMeans(pred))

library(BASS)
mod.normal<-bass(x,y)
plot(mod.normal)
var(f(xtest)-colMeans(predict(mod.normal,xtest)))


dat<-read.csv('~/git/realTime/AlAl/data_trial5/features_cdf104S.csv')
X<-read.table('~/git/realTime/AlAl/Al.trial5.design.txt',head=T)

yy<-dat[,7]
yy<-scale(yy)
mod.normal<-bass(X,yy)
plot(mod.normal)

mod<-tbass(X,yy,nu=10)
plot(sqrt(mod$s2[-c(1:1000)]),type='l')
plot(mod$nbasis,type='l')
hist(colMeans(1/mod$v[-c(1:9000),]))

