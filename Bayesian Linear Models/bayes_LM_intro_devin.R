dat<-read.csv('Desktop/data.csv')
y<-dat$y
X<-as.matrix(cbind(1,dat[,-c(1,2)]))
library(mnormt)

a<-0
b<-0
tau2<-100
t2I<-diag(p)/tau2
XtX<-t(X)%*%X
Xty<-t(X)%*%y

p<-ncol(X)
n<-nrow(X)
nmcmc<-2000
beta<-matrix(nrow=nmcmc,ncol=p)
s2<-rep(NA,nmcmc)
beta[1,]<-0
s2[1]<-1


for(i in 2:nmcmc){
  
  s2[i]<-1/rgamma(1,n/2+a,b+.5*sum((y-X%*%beta[i-1,])^2))
  
  S<-pd.solve(XtX/s2[i]+t2I)
  m<-S%*%Xty/s2[i]
  beta[i,]<-rmnorm(1,m,S)
  
}

burn<-1:1000
matplot(beta,type='l')
plot(sqrt(s2[-burn]),type='l')

colMeans(beta[-burn,])
mean(sqrt(s2[-burn]))
t(apply(beta,2,quantile,probs=c(.025,.975)))


mod<-lm(y~X-1)
mod$coef
summary(mod)$sigma

confint(mod)


