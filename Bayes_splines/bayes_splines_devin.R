set.seed(12)
n<-300
x<-seq(0,1,length.out=n)
y<-sin(2*pi*x^2)*10+rnorm(n)
plot(x,y)

library(mnormt)

a<-0
b<-0
tau2<-100
p<-15
t2I<-diag(p)/tau2

pos<-function(x){
  x[x<0]<-0
  x
}
spline.basis<-function(x,knots,signs,degree=1){
  nb<-length(knots)
  out<-matrix(nrow=length(x),ncol=nb)
  for(i in 1:nb)
    out[,i]<-pos((x-knots[i])*signs[i])^degree
  out
}

nmcmc<-20000
beta<-matrix(nrow=nmcmc,ncol=p)
s2<-rep(NA,nmcmc)
knots<-signs<-matrix(nrow=nmcmc,ncol=p)
beta[1,]<-0
s2[1]<-1
knots[1,]<-runif(p)
signs[1,]<-1
X<-spline.basis(x,knots[1,],signs[1,])
bhat<-pd.solve(crossprod(X)/1+t2I)%*%crossprod(X,y)/1
count=0
for(i in 2:nmcmc){
  
  s2[i]<-1/rgamma(1,n/2+a,b+.5*sum((y-X%*%beta[i-1,])^2))
  
  knots[i,]<-knots[i-1,]
  signs[i,]<-signs[i-1,]
  knots.cand<-runif(p)
  signs.cand<-sample(c(-1,1),size=p,replace=T)
  X.cand<-spline.basis(x,knots.cand,signs.cand)
  bhat.cand<-pd.solve(crossprod(X.cand)/s2[i]+t2I)%*%crossprod(X.cand,y)/s2[i]
  alpha<- -.5/s2[i]*sum((y-X.cand%*%bhat.cand)^2) - (-.5/s2[i]*sum((y-X%*%bhat)^2))
  if(log(runif(1)) < alpha){
    knots[i,]<-knots.cand
    signs[i,]<-signs.cand
    X<-X.cand
    bhat=bhat.cand
    count=count+1
  }
  
  S<-pd.solve(crossprod(X)/s2[i]+t2I)
  m<-S%*%crossprod(X,y)/s2[i]
  beta[i,]<-rmnorm(1,m,S)
  
}

burn<-1:1000
matplot(beta,type='l')
plot(sqrt(s2[-c(1:10)]),type='l')

pred<-matrix(nrow=nmcmc-length(burn),ncol=n)
k=0
for(i in (length(burn)+1):nmcmc){
  k=k+1
  pred[k,]<-spline.basis(x,knots[i,],signs[i,])%*%beta[i,]
}

matplot(x,t(pred),type='l')
points(x,y)
curve(sin(2*pi*x^2)*10,add=T,lwd=3)
