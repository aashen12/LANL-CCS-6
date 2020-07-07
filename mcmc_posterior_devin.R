n<-100
data <- rmvnorm(n, c(1, 2, 3), cbind(c(1, 1.4, 2.1), c(1.4, 4.0, 4.2), c(2.1, 4.2, 9.0)))
p<-3

library(mnormt)
library(invgamma)

s2.a<-.01
s2.b<-.01
rho.a<-1
rho.b<-1
m<-rep(0,3)
Sinv<-diag(3)/100

nybar<-colMeans(data)*n
nmcmc<-3000
mu<-s2<-rho<-matrix(nrow=nmcmc,ncol=p)
mu[1,]<-0
s2[1,]<-.1
rho[1,]<-0

makeCov<-function(cors,vars){
  cor.mat<-diag(length(vars))
  cor.mat[upper.tri(cor.mat)]<-cor.mat[lower.tri(cor.mat)]<-cors
  mat<-matrix(rep(sqrt(vars),length(vars)),nrow=length(vars))
  out<-cor.mat*mat*t(mat)
  out[upper.tri(out)]<-out[lower.tri(out)]
  out
}

lpost<-function(mu.curr,s2.curr,rho.curr){
  Sig<-makeCov(rho.curr,s2.curr)
  if(any(eigen(Sig)$values<0))
    return('fail')
  temp <- t(t(data) - mu.curr) 
  llike<- -.5*sum((temp %*% pd.solve(Sig)) * temp) -n/2*determinant(Sig)$mod
  lprior<- sum(dinvgamma(s2.curr,s2.a,s2.b,log=T)) + sum(dbeta((rho.curr+1)*.5,rho.a,rho.b))
  llike + lprior
}


sd.s2<-c(.5,1.5,2)
sd.rho<-c(.15,.15,.15)
count.s2<-count.rho<-rep(0,p)
for(i in 2:nmcmc){
  #rho[i,]<-.7
  Sig<-makeCov(rho[i-1,],s2[i-1,])
  
  S.star<-pd.solve(Sinv + n*pd.solve(Sig))
  m.star<-S.star%*%(Sinv%*%m + pd.solve(Sig)%*%nybar)
  mu[i,]<-rmnorm(1,m.star,S.star)
  
  s2[i,]<-s2[i-1,]
  for(j in 1:3){
    s2j.cand<-rtruncnorm(1,a=0,mean=s2[i-1,j],sd=sd.s2[j])
    s2.cand<-s2[i,]
    s2.cand[j]<-s2j.cand
    
    if(lpost(mu[i,],s2.cand,rho[i-1,])=='fail'){
      alpha<- -1e10
    } else{
      alpha<-(lpost(mu[i,],s2.cand,rho[i-1,]) - log(dtruncnorm(s2j.cand,a=0,mean = s2[i-1,j],sd = sd.s2[j]))) - 
        (lpost(mu[i,],s2[i,],rho[i-1,]) - log(dtruncnorm(s2[i-1,j],a=0,mean = s2j.cand,sd = sd.s2[j])))
    }
    if(log(runif(1))<alpha){
      s2[i,j]<-s2j.cand
      count.s2[j]<-count.s2[j]+1
    }
  }
  
  
  rho[i,]<-rho[i-1,]
  for(j in 1:3){
    rhoj.cand<-rtruncnorm(1,a=-.99,b=.99,mean=rho[i-1,j],sd=sd.rho[j])
    rho.cand<-rho[i,]
    rho.cand[j]<-rhoj.cand
    
    if(lpost(mu[i,],s2[i,],rho.cand)=='fail'){
      alpha<- -1e10
    } else{
      alpha<-(lpost(mu[i,],s2[i,],rho.cand) - log(dtruncnorm(rhoj.cand,a=-1,b=1,mean = rho[i-1,j],sd = sd.rho[j]))) - 
        (lpost(mu[i,],s2[i,],rho[i,]) - log(dtruncnorm(rho[i-1,j],a=-1,b=1,mean = rhoj.cand,sd = sd.rho[j])))
    }
    if(log(runif(1))<alpha){
      rho[i,j]<-rhoj.cand
      count.rho[j]<-count.rho[j]+1
    }
  }
  
  
  
}

matplot(mu,type='l')
abline(h=1:3)
matplot(sqrt(s2),type='l')
abline(h=1:3)
matplot(rho,type='l')
abline(h=.7)
count.s2/nmcmc
count.rho/nmcmc

burn<-1:1000
colMeans(mu[-burn,])
colMeans(s2[-burn,])
colMeans(rho[-burn,])
