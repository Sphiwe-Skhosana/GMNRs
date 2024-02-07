## The ```locpol``` package is
if(!require(locpol)){install.packages("locpol")}else{library(locpol)}
library(splines)


## A function that computes the BIC for the fitted model
BIC=function(t,x,bw,K,LogLik){
  n=length(t)
  p=ncol(x)
  rk=2.5375;ck=0.7737;
  dfm=(rk*abs(diff(range(t)))*ck)/bw
  dfn=(K-1)*dfm+(K+p*K)
  BIC=-2*LogLik+log(n)*dfn
  return(BIC)
}


##Finding the trimmed (or inter-percentile) range
trim=function(x){
  n=length(x)
  idx=1:(n/2)
  y1=x[idx];y2=sort(x[-idx])
  return(abs(y1-y2))
}

## The kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  f=ifelse(abs(z)<=1, 0.75*(1 - z^2),0)/h
  #f=dnorm(z)
  out=f
  if(sum(f)>0){out=f/sum(f)};
  return(out)
}

## A function that computes the density of a Gaussian mixture model (GMM)
GMM=function(y,mix.mu,mix.prop,mix.sigma){
  k=length(mix.mu)
  out=rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[j],0,mix.sigma[j])))
  return(out)
}

##Local polynomial smoother
local.polynomial.smoother=function(x,y,xgrid,bw,d,W){
  n=length(y)
  g=locPolSmootherC(x,y,xgrid,bw,d,gaussK,weig=W)
  return(g)
}

##Backfitting function
backfit=function(y,t,x,tgrid,k,mh,Beta_init,pi_init,sigma2_init,bw,backfit=TRUE){
  n=length(y)
  z=as.matrix(x)
  u=t
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  Beta0=Beta_init
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  Kh=sapply(tgrid,function(t0) Kern(u,t0,bw))
  diff=1e10
  count=0
  while(diff>1e-10){
    #E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j]-z%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-100)
    gn=g/rowSums(g)
    #M-Step
    pi1=colSums(gn)/n
    sigma21=Beta1=NULL
    for(j in 1:k){
      W=gn[,j]*Kh
      S=diag(W);I=diag(n);R=diag(gn[,j])
      #yt=(I-S)%*%y;zt=(I-S)%*%z
      Beta1=cbind(Beta1,solve(t(z)%*%R%*%z)%*%t(z)%*%R%*%(y-mh[,j]))
      res2=c((y-mh[,j]-z%*%Beta1[,j])^2)
      sigma21=c(sigma21,sum(gn[,j]*res2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j]-z%*%Beta1[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    sigma20=sigma21
    pi0=pi1
    Beta0=Beta1
    count=count+1
    if(count==1e2) diff=1e-100
  }
  
  #Re-estimating the non-parametric functions given the global parameter estimates
  mu=mh
  mu0=mh
  if(backfit){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(tgrid,function(x0) Kern(x,x0,bw))
    ngrid=length(tgrid)
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      #g=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[t,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j]-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      #mu1=t(sapply(1:ngrid,function(t){
      mu1=sapply(1:k,function(j){
        W=gn[,j]
        res=c(y-z%*%Beta1[,j])
        local.polynomial.smoother(u,res,tgrid,bw,0,W)[,2]
      })
      #}))
      mu=sapply(1:k,function(j) approx(tgrid,mu1[,j],u,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j]-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mu
      count=count+1
      if(count==1e2) diff=1e-100
    }
  }
  out=list(mu=mu,pi1=pi1,sigma21=sigma21,Beta1=Beta1)
}

##A function that performs iterative backfitting
backfit_fullyIter=function(y,t,x,tgrid,k,mh,Beta_init,pi_init,sigma2_init,bw){
  n=length(y)
  z=as.matrix(x)
  u=t
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  Beta0=Beta_init
  LogLikn=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  difff=1e10
  countN=0
  while(difff>1e-10){
    LogLik0=LogLikn
    diff=1e10
    count=0
    Kh=sapply(tgrid,function(t0) Kern(u,t0,bw))
    while(diff>1e-10){
      #E-Step
      g=sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j]-z%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-100)
      gn=g/rowSums(g)
      #M-Step
      pi1=colSums(gn)/n
      sigma21=Beta1=NULL
      for(j in 1:k){
        W=gn[,j]*Kh
        S=diag(W);I=diag(n);R=diag(gn[,j])
        #yt=(I-S)%*%y;zt=(I-S)%*%z
        Beta1=cbind(Beta1,solve(t(z)%*%R%*%z)%*%t(z)%*%R%*%(y-mh[,j]))
        res2=c((y-mh[,j]-z%*%Beta1[,j])^2)
        sigma21=c(sigma21,sum(gn[,j]*res2)/sum(gn[,j]))
      }
      #Evaluate for convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j]-z%*%Beta1[,j],0,sqrt(sigma21[j]))))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      sigma20=sigma21
      pi0=pi1
      Beta0=Beta1
      count=count+1
      if(count==5e2) diff=1e-100
    }
    #Re-estimating the non-parametric functions given the global parameter estimates
    mu0=mh
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j]-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(tgrid,function(u0) Kern(u,u0,bw))
    ngrid=length(tgrid)
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      #g=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[t,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j]-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      #mu1=t(sapply(1:ngrid,function(t){
      mu1=sapply(1:k,function(j){
        W=gn[,j]
        res=c(y-z%*%Beta1[,j])
        local.polynomial.smoother(u,res,tgrid,bw,0,W)[,2]
      })
      #}))
      mu=sapply(1:k,function(j) approx(tgrid,mu1[,j],u,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j]-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mu1
      count=count+1
      if(count==5e2) diff=1e-100
    }
    ##Evaluate for overall convergence
    difff=abs(LogLik1-LogLikn)
    LogLikn=LogLik1
    mh=mu1
    Beta0=Beta1
    pi0=pi1
    sigma20=sigma21
    countN=countN+1
    if(countN==1e2) difff=1e-100
  }
  out=list(mu=mu,pi1=pi1,sigma21=sigma21,Beta1=Beta1)
}

##A function used to initialize the algorithms that fit the SPGMPLMs
initialize.model=function(x,t,y,k,method=NULL,true.init=NULL,p){
  n=length(y)
  BIC=1e10
  if(method==1){##Mixture of partial regression splines
    for(j in 1:1e2){
      m=list(BIC=1e6)
      try({m=mix.partial.reg.splines(x,t,y,k)},silent=T)
      if(m$BIC<BIC){init.model=m$init.model0;BIC=m$BIC}
    }
  }
  if(method==2){##Mixture of linear regression
    for(j in 1:1e2){
      m=list(BIC=1e6)
      try({m=mix.poly(x,y,2,p)})
      if(m$BIC<BIC){init.model=m$init.model;BIC=m$BIC}
    }
  }
  if(method==3){##True values
    m0=true.init
    init.model=list(mu0=m0$mu,Beta0=t(as.matrix(m0$beta)),sigma20=m0$sigma2,pi0=m0$rho)
  }
  return(list(init.model0=init.model))
}

##A function that computes the roughness of a curve
Rough_curve<-function(x,f=NULL){
  if(!is.null(f)) f=f[order(x)]
  x=sort(x)
  n=length(x);
  N=diag(n)
  ##Checking whether the x's are distinct
  if(length(unique(x))!=n){z=unique(x);n=length(z);N=sapply(z,function(u) as.numeric(x==u));  x=z}
  h<-NULL
  for(i in (1:(n-1))){h=c(h,x[i+1]-x[i])}
  Q<-matrix(0,nrow=n,ncol=n-2)
  for(j in (2:(n-1))){Q[j-1,j-1]=1/h[j-1];Q[j,j-1]=-(1/h[j-1]+1/h[j]);Q[j+1,j-1]=1/h[j]}
  R<-matrix(0,n-2,n-2)
  for(i in (2:(n-1))){R[i-1,i-1]=(h[i-1]+h[i])/3;if(i<=(n-2)){R[i,i-1]=R[i-1,i]=h[i]/6}}
  K=Q%*%solve(R)%*%t(Q)
  Rh=NULL
  if(!is.null(f)){
    f=apply(N,2,function(x) mean(f[which(x==1)])); 
    Rh=t(f)%*%K%*%f} ##Roughness value
  return(list(Rh=Rh,K=K,N=N))
}

###A function that fits the mixture of partial regression splines (to initialize) 
mix.partial.reg.splines=function(X,t,y,k){
  library(splines)
  X=as.matrix(X)
  n=nrow(X)
  Xd=X##Design matrix
  Td=bs(t,knots=quantile(t,probs=seq(0.25,0.75,0.25)),intercept=T)##B-spline basis matrix
  p=ncol(X);d=ncol(Td);
  ##Initial state
  pi0=rep(1/k,k)
  sigma20=rgamma(k,1,1)^2
  Beta0=matrix(rnorm(k*p),p,k)
  lmd0=matrix(rnorm(k*d),d,k)
  LogLik0=sum(log(rowSums(pi0*dnorm(y-Xd%*%Beta0-Td%*%lmd0,0,sqrt(sigma20)))))
  diff=1e6
  ll=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    g=pi0*dnorm(y-Xd%*%Beta0-Td%*%lmd0,0,sqrt(sigma20))
    gn=g/rowSums(g)
    ##M-Step
    pi1=colSums(gn)/n
    lmd1=Beta1=sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      res0=as.numeric(y-Td%*%lmd0[,j])
      Beta1=cbind(Beta1,solve(t(Xd)%*%W%*%Xd)%*%t(Xd)%*%W%*%res0)
      res1=as.numeric(y-Xd%*%Beta1[,j])
      lmd1=cbind(lmd1,solve(t(Td)%*%W%*%Td)%*%t(Td)%*%W%*%res1)
      sigma21=c(sigma21,sum(gn[,j]*(y-Xd%*%Beta1[,j]-Td%*%lmd1[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(pi1*dnorm(y-Xd%*%Beta1-Td%*%lmd1,0,sqrt(sigma21)))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    Beta0=Beta1
    lmd0=lmd1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==1e3) diff=1e-100
  }
  df=k*(p+d)+2*k-1
  df_reg=k*(p+d)
  BIC=-2*LogLik0+df*log(n)
  mu1=Td%*%lmd1
  r0=pi1*dnorm(y-Xd%*%Beta1-Td%*%lmd1,0,sqrt(sigma21));r=r0/rowSums(r0)
  model0=list(r=r,Beta0=matrix(Beta1,p,k),pi0=pi1,sigma20=sigma21,mu0=mu1)
  return(list(r=r,BIC=BIC,lmd0=lmd1,mu=mu1,Beta=Beta1,pi=pi1,sigma2=sigma21,init.model0=model0,df_reg=df_reg))
}

## A function that fits a Gaussian mixture of linear regressions (GMLRs) model
GaussLinMix=function(x,y,x0,k,weights=NULL,model0){
  n=length(y)
  x=as.matrix(x)
  z=cbind(rep(1,n),x)
  pi0=model0$pi0
  sigma20=model0$sigma20
  Beta0=rbind(model0$mu0,model0$Beta0)
  LL0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j]))+1e-100;gn=g/rowSums(g)
    ##M-step
    pi1=colSums(weights*gn)/sum(weights)
    Beta1=sigma21=mu1=NULL
    for(j in 1:k){
      W=weights*gn[,j]
      Beta=solve(t(z)%*%diag(W)%*%z)%*%t(z)%*%diag(W)%*%y
      Beta1=cbind(Beta1,Beta)
      sigma21=c(sigma21,sqrt(sum(W*(y-z%*%Beta)^2)/sum(W)))
    }
    
    LL1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
    diff=abs(LL1-LL0)
    LL0=LL1
    Beta0=Beta1
    pi0=pi1
    count=count+1
    sigma20=sigma21
    llk=c(llk,LL0)
    if(count==1e3) diff=1e-100
  }
  mu1=Beta1[1,];Beta1=Beta1[-1,]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  return(list(resp=gn,mix.mu=mu1,mix.Beta=Beta1,mix.sigma2=sigma21,mix.prop=pi1,LL=LL1))
}

##A function that fits the SPGMPLMs using the profile-likelihood EM (PL-EM) algorithm
SPGMPLMs_PL_EM=function(x,t,y,k,bw,tgrid,init.model){
  n=length(y)
  z=as.matrix(x)
  ngrid=length(tgrid)
  ##Initial state
  pi0=init.model$pi0
  mu0=init.model$mu0;
  Beta0=init.model$Beta0
  sigma20=init.model$sigma20
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  Kh=sapply(tgrid,function(t0) Kern(t,t0,bw))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])+1e-100);gn=g/rowSums(g)
    ##M-step
    pi1=colSums(gn)/n
    Beta1=mu1=sigma21=NULL
    for(j in 1:k){
      W=gn[,j];R=diag(gn[,j])
      res=c(y-z%*%Beta0[,j])
      mh=local.polynomial.smoother(t,res,tgrid,bw,0,W)[,2]
      mu1=cbind(mu1,approx(tgrid,mh,t,rule=2)$y)
      #S=diag(W);I=diag(n)
      #xt=(I-S)%*%z;yt=(I-S)%*%y
      Beta1=cbind(Beta1,solve(t(z)%*%R%*%z)%*%t(z)%*%R%*%(y-mu1[,j]))
      sigma21=cbind(sigma21,sum(gn[,j]*(y-mu1[,j]-z%*%Beta1[,j])^2)/sum(gn[,j]))
    }
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j]-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    Beta0=Beta1
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  res=order(Beta1[1,]);mu1=mu1[,res];Beta1=t(as.matrix(Beta1[,res]));pi1=pi1[res];sigma21=sigma21[res]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j]-z%*%Beta1[,j],0,sqrt(sigma21[j])));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.beta=Beta1,mix.sigma2=sigma21,df=df,LL=LL1)
  return(out)
}

##A function that fits the SPGMPLMs using the proposed objective-based profile-likelihood EM (OB-PL-EM) algorithm
SPGMPLMs_OB_PL_EM=function(x,t,y,k,bw,tgrid,init.model){
  n=length(y)
  z=as.matrix(x);p=ncol(z)
  u=t
  ngrid=length(tgrid)
  ##Initial State
  pi0=init.model$pi0
  mu0=init.model$mu0;
  Beta0=init.model$Beta0;
  sigma20=init.model$sigma20
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mu0[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],tgrid,rule=2)$y)
  Kh=sapply(tgrid,function(u0) Kern(u,u0,bw))
  diff=1e6
  count=0
  R0=1e5
  llk=NULL
  while(diff>1e-10){
    ##E-Step
    zh=lapply(1:ngrid,function(t){
      model0=list(pi0=pi0,sigma20=sigma20,mu0=mu0[t,],Beta0=Beta0)
      model=GaussLinMix(z,y,tgrid[t],k,Kh[,t],model0);list(resp=model$resp,LL=model$LL)})
    
    ##M-Step
    mu=lapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        W=zh[[t]]$resp[,j]
        yh=z%*%Beta0[,j] ##linear (parametric) term
        res=c(y-yh)
        mh=local.polynomial.smoother(u,res,tgrid,bw,0,W)[,2]
 ##non-parametric term
      })
    })
    ##Measuring the roughness
    #plot(u,y)
    for(t in 1:ngrid){
      #for(j in 1:k) lines(tgrid,mu[[t]][,j],col="red")
      R00=NULL
      for(j in 1:k) {R00=c(R00,Rough_curve(tgrid,mu[[t]][,j])$Rh)};R1=sum(log(R00))
      dist=max(apply(round(mu[[t]],4),1,function(x) abs(diff(x))))
      D=apply(mu[[t]],2,function(x) abs(trim(x)))
      id=apply(D,2,sum);id=which.min(id)
      ##Model space
      #No CRF can be a straight line (D>>0); two CRFs cannot be equal(dist>>0) 
      #How can I restrict the admission of straight lines? (sol:not all trimmed ranges should be zero)
      if(R1<R0 & dist>1e-1 & !all(D<0.5) & !all(D[,id]==0)){
        #if(R1<R0){
        mu1=mu[[t]];
      }
    }
    ##Evaluating convergence
    mix.mu=NULL
    for(j in 1:k){
      mix.mu=cbind(mix.mu,approx(tgrid,mu1[,j],xout=u,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mix.mu[,j]-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
    diff=abs(LogLik1-LogLik0)
    mu0=mu1
    LogLik0=LogLik1
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,u,x,tgrid,k,mu1,Beta0,pi0,sigma20,bw);
  mu2=out$mu;Beta1=out$Beta1;pi1=out$pi1;sigma21=out$sigma21;
  res=order(Beta1[1,]);mu2=mu2[,res];Beta1=t(as.matrix(Beta1[,res]));pi1=pi1[res];sigma21=sigma21[res]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j]-z%*%Beta1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.beta=Beta1,mix.sigma2=sigma21,LL=LL1,R=R1)
  return(out)
}

## A function that fits the SPGMPLMs using the Naive EM algorithm
SPGMPLMs_Naive_EM=function(X,t,y,k,bw,tgrid,init.model=NULL){
  X=as.matrix(X)
  u=t
  p=ncol(X);n=nrow(X)
  ngrid=length(tgrid)
  ##Initial state
  pi0=init.model$pi0
  Beta0=init.model$Beta0; if(!is.matrix(Beta0)) Beta0=t(as.matrix(Beta0))
  mu0=init.model$mu0
  sigma20=init.model$sigma20
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-X%*%Beta0[,j]-mu0[,j],0,sqrt(sigma20)[j])))))
  Kh=sapply(tgrid,function(u0) Kern(u,u0,bw))
  diff=1e6
  llk=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    g=lapply(1:ngrid,function(i){g=sapply(1:k,function(j) pi0[j]*dnorm(y-X%*%Beta0[,j]-mu0[i,j],0,sqrt(sigma20[j])));gn=g/rowSums(g)})
    ##M-Step
    mu=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        W=g[[t]][,j]
        yh=X%*%Beta0[,j] ##linear (parametric) term
        res=c(y-yh)
        mh=local.polynomial.smoother(u,res,tgrid[t],bw,0,W)[,2]

      })
    }))
    ##Evaluating convergence
    mix.mu=NULL
    for(j in 1:k){
      mix.mu=cbind(mix.mu,approx(tgrid,mu[,j],xout=u,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mix.mu[,j]-X%*%Beta0[,j],0,sqrt(sigma20)[j])))))
    diff=abs(LogLik1-LogLik0)
    mu0=mu
    LogLik0=LogLik1
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,u,X,tgrid,k,mu1,Beta0,pi0,sigma20,bw,backfit=F);
  mu2=out$mu;Beta1=out$Beta1;pi1=out$pi1;sigma21=out$sigma21;
  res=order(Beta1[1,]);mu2=mu2[,res];Beta1=t(as.matrix(Beta1[,res]));pi1=pi1[res];sigma21=sigma21[res]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j]-X%*%Beta1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.mu=mu2,mix.prop=pi1,mix.beta=Beta1,mix.sigma2=sigma21,LL=LL1)
  return(out)
}
