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
k3=step=1

i4= 178
i5=    146
thre1=   1.299e+03
tthre1= 3.388e+03
k=17
#############################################
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
################################################################################ calculate P-value for C0
fit1=lm(crcase~lag1crcase+sum+r+thre3+t+tthre3+0) 
fit2=lm(crcase~           sum+r+thre3+t+tthre3+0)
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
###################################################################################################################
X=cbind(summ,r1,th1,t1,tth1     )
beta=fit2$coefficients
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
      X[,1][(k+i+1):length(X[,2])]=X[,1][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lag1crcase[(k+i+1):(length(y_hat) )]
    } else if ( i>k & i!=((length(crcase)-k  ) )) {
      
      lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
      X[,1][(k+i+1):length(X[,2])]=X[,1][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lag1crcase[(k+i+1):(length(y_hat) )]
    }   
  }
  return(err/(length(crcase)-k) *((length(crcase)-k) +1)/2)
}
est.b1 = optim(beta, obj,method="CG")

beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[2]*(r1)[i]*y[i-1]+
        beta[3]*th1[i]*y[i-1]+
        beta[4]*(t1)[i]*y[i-1]+
        beta[5]*(tth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b11 = optim(beta, obj)




beta=est.b$par
errz=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*th1[i]*y[i-1]+
      beta[5]*(t1)[i]*y[i-1]+
      beta[6]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  errz = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+errz
}
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[2]*(r1)[i]*y[i-1]+
      beta[3]*th1[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+beta[5]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}


betax=est.b11$par
T07=(err1-errz)/errz

cl = makeCluster(num_cores)
registerDoParallel(cl)
KK=1000
########################## find the best parameters and you can choose your own K
ress7 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[2]*(r1)[i]*y[i-1]+
      beta[3]*th1[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+beta[5]*(tth1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~      sum+r+thre3+t+tthre3+0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #########################################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  #######################################
  X=cbind(summ,br1,bth1,bt1,btth1     )
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,1][(k+i+1):length(X[,2])]=X[,1][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,1][(k+i+1):length(X[,2])]=X[,1][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  
  
  
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[2]*(br1)[i]*y[i-1]+
              beta[3]*bth1[i]*y[i-1]+
              beta[4]*(bt1)[i]*y[i-1]+
              beta[5]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress7
TB7=res[,1]
################################################################################# calculate P-value for CI
fit2=lm(crcase~lag1crcase+r+thre3+t+tthre3+0)
#################################################################################
X=cbind(rep(1,length(crcase)),r1,th1,t1,tth1     )
beta=fit2$coefficients
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
      #X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lag1crcase[(k+i+1):(length(y_hat) )]
    } else if ( i>k & i!=((length(crcase)-k  ) )) {
      
      lag1crcase[(k+i+1):length(lag1crcase)]=y_hat[(k+i):(length(y_hat)-1)]
      #X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lag1crcase[(k+i+1):(length(y_hat) )]
    }   
  }
  return(err/(length(crcase)-k) *((length(crcase)-k) +1)/2)
}
est.b1 = optim(beta, obj,method="CG")
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+
      beta[2]*(r1)[i]*y[i-1]+
      beta[3]*(th1)[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}

beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*y[i-1]+
        beta[2]*(r1)[i]*y[i-1]+
        beta[3]*th1[i]*y[i-1]+
        beta[4]*(t1)[i]*y[i-1]+
        beta[5]*(tth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b12 = optim(beta, obj)


betax=est.b12$par
T08=(err1-errz)/errz



cl = makeCluster(num_cores)
registerDoParallel(cl)
ress8 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*y[i-1]+
      beta[2]*(r1)[i]*y[i-1]+
      beta[3]*th1[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+beta[5]*(tth1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~lagy1+    r+thre3+t+tthre3+0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #######################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  
  X=cbind(rep(1,length(y1)),br1,bth1,bt1,btth1     )
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        #X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        #X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  
  
  
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=   beta[1]*y[i-1]+
        beta[2]*(br1)[i]*y[i-1]+
        beta[3]*bth1[i]*y[i-1]+
        beta[4]*(bt1)[i]*y[i-1]+
        beta[5]*(btth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress8
TB8=res[,1]


################################################################################# calculate P-value for Cp
fit2=lm(crcase~lag1crcase+sum+thre3+t+tthre3+0)
#################################################################
X=cbind(rep(1,length(crcase)),summ,th1,t1,tth1     )
beta=fit2$coefficients
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
est.b1 = optim(beta, obj,method="CG")
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(th1)[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}


beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*th1[i]*y[i-1]+
        beta[4]*(t1)[i]*y[i-1]+
        beta[5]*(tth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b13 = optim(beta, obj)


betax=est.b13$par
T09=(err1-errz)/errz
cl = makeCluster(num_cores)
registerDoParallel(cl)
ress9 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(th1)[i]*y[i-1]+
      beta[4]*t1[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~lagy1+sum+  thre3+t+tthre3+0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #######################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  
  X=cbind(rep(1,length(y1)),summ,bth1,bt1,btth1     )
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  
  
  
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*bth1[i]*y[i-1]+
              beta[4]*(bt1)[i]*y[i-1]+
              beta[5]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress9
TB9=res[,1]

################################################################################# calculate P-value for Cp'
fit2=lm(crcase~lag1crcase+sum+r+t+tthre3+0)
###################################################################################################################
X=cbind(rep(1,length(crcase)),summ,r1,t1,tth1     )
beta=fit2$coefficients
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
est.b1 = optim(beta, obj,method="CG")
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*(t1)[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}


beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*(r1)[i]*y[i-1]+
        beta[4]*(t1)[i]*y[i-1]+
        beta[5]*(tth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b14 = optim(beta, obj)
betax=est.b14$par
T010=(err1-errz)/errz



cl = makeCluster(num_cores)
registerDoParallel(cl)
ress10 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*t1[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~lagy1+sum+r      +t+tthre3+0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #######################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  
  X=cbind(rep(1,length(y1)),summ,br1,bt1,btth1     )
  
  
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  
  
  
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*br1[i]*y[i-1]+
        beta[4]*(bt1)[i]*y[i-1]+
        beta[5]*(btth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress10
TB10=res[,1]
################################################################################# calculate P-value for Ct
fit2=lm(crcase~lag1crcase+sum+r+thre3+tthre3+0)
###################################################################################################################
X=cbind(rep(1,length(crcase)),summ,r1,th1,tth1     )
beta=fit2$coefficients
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
est.b1 = optim(beta, obj,method="CG")
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*(th1)[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}


beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*(r1)[i]*y[i-1]+
        beta[4]*(th1)[i]*y[i-1]+
        beta[5]*(tth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b15 = optim(beta, obj)


betax=est.b15$par
T011=(err1-errz)/errz


cl = makeCluster(num_cores)
registerDoParallel(cl)
ress11 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*th1[i]*y[i-1]+
      beta[5]*(tth1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~lagy1+sum+r+thre3  +tthre3+0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #######################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  
  X=cbind(rep(1,length(y1)),summ,br1,bth1,btth1     )
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  
  
  
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*br1[i]*y[i-1]+
        beta[4]*(bth1)[i]*y[i-1]+
        beta[5]*(btth1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress11
TB11=res[,1]
################################################################################# calculate P-value for Ct'
fit2=lm(crcase~lag1crcase+sum+r+thre3+t+0)
###############################################################################
X=cbind(rep(1,length(crcase)),summ,r1,th1,t1     )
beta=fit2$coefficients
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
est.b1 = optim(beta, obj,method="CG")
################################################################
beta=est.b1$par
err1=0
for (k1 in  seq( (k),(length(crcase)-k3),1 )    ) {
  y=crcase[1:k1]
  for(i in (k1+1):(length(crcase)))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*(th1)[i]*y[i-1]+
      beta[5]*(t1)[i]*y[i-1]
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  err1 = mean((crcase[(k1+1):(length(crcase))]-y[(k1+1):(length(crcase))]   )^2)+err1
}


beta=fit2$coefficients 
obj = function(beta)
{
  err=0
  for (k1 in  seq( (k),(length(crcase)-k3),step )    ) {
    y=crcase[1:k1]
    for(i in (k1+1):(k1+k3))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*(r1)[i]*y[i-1]+
        beta[4]*(th1)[i]*y[i-1]+
        beta[5]*(t1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err = -sum(  crcase[(k1+1):(k1+k3)]*  log(y)[(k1+1):(k1+k3)]   -y[(k1+1):(k1+k3)]      )+err
  }
  err=err/length(     seq( (k),(length(crcase)-k3),step )            )
  return(err)
}
est.b16 = optim(beta, obj)
betax=est.b16$par


T012=(err1-errz)/errz

cl = makeCluster(num_cores)
registerDoParallel(cl)
ress12 = foreach(iter = 1:KK,.combine='rbind')%dopar%{
  beta=betax
  y=crcase[1:k]
  for(i in (k+1):length(crcase))
  {
    y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
      beta[3]*(r1)[i]*y[i-1]+
      beta[4]*th1[i]*y[i-1]+
      beta[5]*(t1)[i]*y[i-1]
    y[i]=rpois(1,y[i])
    y[i]=max(1,y[i])
    y[i]=min(2*max(crcase),y[i])
  }
  
  
  y1=y[(k+1):length(crcase)]
  lagy1=y[(k):(length(crcase)-1)]
  lagy30=c(  rep(1,30-k), y     )
  summ=numeric()
  for(i in 1:length(y1))
  {summ[i]=sum(lagy30[(30+i-k):(29+i)])}
  sum=summ*lagy1
  
  r= (rain_1[i4,][(k+1):length(rain_1[i4,])])*lagy1
  t= (temp_1[i5,][(k+1):length(rain_1[i4,])])*lagy1
  
  
  
  
  
  thre2=((rain_1[i4,][(k+1):length(rain_1[i4,])])-thre1)
  thre3=thre2*(thre2>0)*lagy1
  
  tthre2=((temp_1[i5,][(k+1):length(rain_1[i4,])])-tthre1)
  tthre3=tthre2*(tthre2>0)*lagy1
  
  fit1=lm(y1~lagy1+sum+r+thre3+t+tthre3+0) ################ result based on 1-step prediction
  fit2=lm(y1~lagy1+sum+r+thre3+t       +0)
  
  br1=(rain_1[i4,][(k+1):length(rain_1[i4,])])
  bt1=(temp_1[i5,][(k+1):length(rain_1[i4,])])
  bth1=thre2*(thre2>0)
  btth1=tthre2*(tthre2>0)
  
  
  
  #######################################################################################
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1,btth1     )
  beta=fit1$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b = optim(beta, obj,method="CG")
  
  
  
  
  ####################################################################################### 
  X=cbind(rep(1,length(y1)),summ,br1,bth1,bt1     )
  beta=fit2$coefficients
  Z=lagy1
  obj = function(beta)
  {
    err=0
    p_y_hat=list()
    for (i in 1:(length(y1)-k)     ) {
      y_hat=(X*lagy1) %*%beta
      y_hat[y_hat<=1]=1
      y_hat[y_hat>=(2*max(y1))]=2*max(y1)
      err=-sum( ( y1[(k+i):length(y1)]*log(y_hat)[(k+i):length(y1)]-y_hat[(k+i):length(y1)] ) )+err
      p_y_hat[[i]]=y_hat
      
      if (i<=k ) {
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-Z[(i+1):(length(Z)-(k) ) ]+lagy1[(k+i+1):(length(y_hat) )]
      } else if ( i>k & i!=((length(y1)-k  ) )) {
        
        lagy1[(k+i+1):length(lagy1)]=y_hat[(k+i):(length(y_hat)-1)]
        X[,2][(k+i+1):length(X[,2])]=X[,2][(k+i):(length(X[,2])-1)]-p_y_hat[[i-k]][i:(length(p_y_hat[[1]])-k-1)]+lagy1[(k+i+1):(length(y_hat) )]
      }   
    }
    return(err/(length(y1)-k) *((length(y1)-k) +1)/2)
  }
  est.b1 = optim(beta, obj,method="CG")
  ####################################################
  beta=est.b$par
  err=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=(beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
              beta[3]*(br1)[i]*y[i-1]+
              beta[4]*bth1[i]*y[i-1]+
              beta[5]*(bt1)[i]*y[i-1]+
              beta[6]*(btth1)[i]*y[i-1])
      y[i]=max(1,y[i])
      y[i]=min(2*max(y1),y[i])
    }
    err = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err
  }
  ################################################################
  beta=est.b1$par
  err1=0
  for (k1 in  seq( (k),(length(y1)-k3),1 )    ) {
    y=y1[1:k1]
    for(i in (k1+1):(length(y1)))
    {
      y[i]=beta[1]*y[i-1]+beta[2]*sum(y[(i-k):(i-1)])*y[i-1]+
        beta[3]*(br1)[i]*y[i-1]+
        beta[4]*bth1[i]*y[i-1]+
        beta[5]*(bt1)[i]*y[i-1]
      y[i]=max(1,y[i])
      y[i]=min(2*max(crcase),y[i])
    }
    err1 = mean((y1[(k1+1):(length(y1))]-y[(k1+1):(length(y1))]   )^2)+err1
  }
  
  T00=(err1-err)/err
  c(T00)
  
}
stopCluster(cl)
res=ress12
TB12=res[,1]
csres1=c( (sum(T07<ress7[,1]))/KK ,(sum(T08<ress8[,1]))/KK ,(sum(T09<ress9[,1]))/KK ,(sum(T010<ress10[,1]))/KK ,
          (sum(T011<ress11[,1]))/KK ,(sum(T012<ress12[,1]))/KK)