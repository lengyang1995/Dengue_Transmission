library(stats)
library(lubridate)
library(scales)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(tidyverse)
########################## extract dengue case count
########################## Please add all files into your R, using Costa Rica as an example.
crcase=SouthAmerica_case$Costa_rica[31:365]
lag1crcase=SouthAmerica_case$Costa_rica[30:364]
lag30crcase=SouthAmerica_case$Costa_rica[1:364]
######################### adjust under-reporting rate
crcase=crcase*25
lag1crcase=lag1crcase*25
lag30crcase=lag30crcase*25
############# set rain and temperature 
temp=Costarica_weather[185:2705,]$temp
rain=Costarica_weather[185:2705,]$rain

h=cumsum(temp)
h=c(0,h)
h1=cumsum(rain)
h1=c(0,h1)

temp_1=matrix(data = NA,nrow=180,ncol=length(crcase))
rain_1=matrix(data = NA,nrow=180,ncol=length(crcase))
for(j in 1:180)
{ 
  for(j1 in 1:length(crcase)){
    temp_1[j,j1]=h[181+7*(j1-1)]-h[181+7*(j1-1)-j]
    rain_1[j,j1]=h1[181+7*(j1-1)]-h1[181+7*(j1-1)-j]
  }
}
############### find the best parameters, K can be chosen from any value
library(stats)
library(doParallel)
library(parallel)
library(foreach)
num_cores = detectCores()
cl = makeCluster(num_cores)
registerDoParallel(cl)
K=10000
########################## find the best parameters and you can choose your own K
ress1 = foreach(iter = 1:K,.combine='rbind')%dopar%{
  k=sample(c(4:30),1)
  i4= sample(c(14:180),1)
  i5=sample(c(14:180),1)
  
  rangeupper=quantile(rain_1[i4,], probs = seq(0.7, 0.8, 0.1))[1]
  rangelower=quantile(rain_1[i4,], probs = seq(0, 0.1, 0.1))[2]
  thre=seq(rangelower, rangeupper, by=5)
  thre1=sample(thre,1)
  
  rangeupper1=quantile(temp_1[i5,], probs = seq(0.7, 0.8, 0.1))[1]
  rangelower1=quantile(temp_1[i5,], probs = seq(0.1, 0.2, 0.1))[2]
  tthre=seq(rangelower1, rangeupper1, by=2)
  tthre1=sample(tthre,1)
  
  
  
  summ=numeric()
  for(i in 1:length(crcase))
  {summ[i]=sum(lag30crcase[(30+i-k):(29+i)])}
  sum=summ*lag1crcase
  
  
  r= (rain_1[i4,])*lag1crcase
  t= (temp_1[i5,])*lag1crcase
  
  
  thre2=(rain_1[i4,]-thre1)
  thre3=thre2*(thre2>0)*lag1crcase
  
  tthre2=(temp_1[i5,]-tthre1)
  tthre3=tthre2*(tthre2>0)*lag1crcase
  
  
  
  fit1=lm(crcase~lag1crcase+sum+r+thre3+t+tthre3+0)

  r1=(rain_1[i4,])
  t1=(temp_1[i5,])
  th1=thre2*(thre2>0)
  tth1=tthre2*(tthre2>0)
  
  X=cbind(rep(1,length(crcase)),summ,r1,th1,t1,tth1     )
  
  ######################## In here, we apply all-step MLE and 'CG'.
  beta=fit1$coefficients
  Z=lag1crcase
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(crcase)-k)     ) {
      y_hat=(X*lag1crcase) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(crcase))]=2*max(crcase)
      err=-sum( ( crcase[(k+i):length(crcase)]*log(y_hat)[(k+i):length(crcase)]-y_hat[(k+i):length(crcase)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lag1crcase[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(crcase)-k  ) )) {
        
        lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lag1crcase[(k+i+1):(length(y_hat) )]
      }  
    }
    return(err/(length(crcase)-k) *((length(crcase)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  c(est.b$value,i4,i5,thre1,tthre1,k)
}
stopCluster(cl)
ress1[which.min(ress1[,1]),]

######## Those are the best parameters from our analysis
i4= 178
i5=    146
thre1=   1.299e+03
tthre1= 3.388e+03
k=17


summ=numeric()
for(i in 1:length(crcase))
{summ[i]=sum(lag30crcase[(30+i-k):(29+i)])}
sum=summ*lag1crcase


r= (rain_1[i4,])*lag1crcase
t= (temp_1[i5,])*lag1crcase


thre2=(rain_1[i4,]-thre1)
thre3=thre2*(thre2>0)*lag1crcase

tthre2=(temp_1[i5,]-tthre1)
tthre3=tthre2*(tthre2>0)*lag1crcase



fit1=lm(crcase~lag1crcase+sum+r+thre3+t+tthre3+0)

r1=(rain_1[i4,])
t1=(temp_1[i5,])
th1=thre2*(thre2>0)
tth1=tthre2*(tthre2>0)

X=cbind(rep(1,length(crcase)),summ,r1,th1,t1,tth1     )


beta=fit1$coefficients
Z=lag1crcase
obj = function(beta)
{
  err=0
  p_y_hat=list()
  for (i in 1:(length(crcase)-k)     ) {
    y_hat=(X*lag1crcase) %*%beta
    y_hat[y_hat<=1]=1
    y_hat[y_hat>=(2*max(crcase))]=2*max(crcase)
    err=-sum( ( crcase[(k+i):length(crcase)]*log(y_hat)[(k+i):length(crcase)]-y_hat[(k+i):length(crcase)] ) )+err
    p_y_hat[[i]]=y_hat
    if (i<=k ) {
      lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
      X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lag1crcase[(k+i+1):(length(y_hat) )]
    } else if ( i>k & i!=((length(crcase)-k  ) )) {
      lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
      X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lag1crcase[(k+i+1):(length(y_hat) )]
    }  
  }
  return(err/(length(crcase)-k) *((length(crcase)-k) +1)/2)
}
est.b = optim(beta, obj,method="CG")
beta = est.b$par

############ re-estimated parameters based on multi-step prediction
y=crcase[1:k]
for(i in (k+1):length(crcase))
{
  y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
    beta[3]*(r1)[i]*y[i-1]+
    beta[4]*(th1)[i]*y[i-1]+
    beta[5]*(t1)[i]*y[i-1]+
    beta[6]*(tth1)[i]*y[i-1]
  y[i]=max(1,y[i])
  y[i]=min(2*max(crcase),y[i])
}


####################### draw graph and see the fitting, which is pretty good.
plot(crcase)
lines(y,col='red')


