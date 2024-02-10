## A function to compute the mixture density of Gaussian mixture model (GMM)
GMM_density=function(y,mix.mu,mix.prop,mix.sigma){
  k=length(mix.mu)
  out=rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[j],0,mix.sigma[j])))
  return(out)
}

##Finding the trimmed (or inter-percentile) range
trim=function(x){
  n=length(x)
  idx=1:(n/2)
  y1=x[idx];y2=sort(x[-idx])
  return(abs(y1-y2))
}

## A kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  f=ifelse(abs(z)<=1, 0.75*(1 - z^2),0)/h
  #f=dnorm(z)
  out=f
  if(sum(f)>0){out=f/sum(f)};
  return(out)
}

##Refitting function
refit=function(y,x,xgrid,k,init.mu,init.prop,init.sigma2,bw,backfit=TRUE){
  n=length(y)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=init.prop
  sigma20=init.sigma2
  mu0=init.mu
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])+1e-100);gn=g/rowSums(g)
    ##M-step
    mu1=pi1=sigma21=NULL
    for(j in 1:k){
      W=gn[,j]*Kh
      pi1=cbind(pi1,colSums(W)/colSums(Kh))
      mh=colSums(W*y)/colSums(W)
      mu1=cbind(mu1,approx(xgrid,mh,x,rule=2)$y)
      sigma21=cbind(sigma21,colSums(W*(y-mu1[,j])^2)/colSums(W))
    }
    pi1=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  out=list(mu=mu1,pi1=pi1,sigma21=sigma21)
}

##Backfitting function
refit_fullyIter=function(y,x,xgrid,k,mh,pi_init,sigma2_init,bw){
  n=length(y)
  ngrid=length(xgrid)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  LogLikn=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mh[,j],0,sqrt(sigma20)[,j])))))
  difff=1e6
  countN=0
  while(difff>1e-10){
    mu0=mh
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
    Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
    diff=1e6
    count=0
    llk=NULL
    while(diff>1e-10){
      ##E-step
      g=sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])+1e-100);gn=g/rowSums(g)
      ##M-step
      mu1=pi1=sigma21=NULL
      for(j in 1:k){
        W=gn[,j]*Kh
        pi1=cbind(pi1,colSums(W)/colSums(Kh))
        mh=colSums(W*y)/colSums(W)
        mu1=cbind(mu1,approx(xgrid,mh,x,rule=2)$y)
        sigma21=cbind(sigma21,colSums(W*(y-mu1[,j])^2)/colSums(W))
      }
      pi1=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
      sigma21=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mu1
      pi0=pi1
      sigma20=sigma21
      count=count+1
      llk=c(llk,LogLik0)
      if(count==1e2) diff=1e-100
    }
    mu1=sapply(1:k,function(j) approx(xgrid,mu0[,j],x,rule=2)$y)
    ##Evaluate for overall convergence
    difff=abs(LogLik1-LogLikn)
    LogLikn=LogLik1
    mh=mu1
    pi0=pi1
    sigma20=sigma21
    countN=countN+1
    if(countN==1e2) difff=1e-100
  }
  out=list(mu=mu1,pi1=pi1,sigma21=sigma21)
}

###Computing the roughness of a curve
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

###A function to initialize the fitting algorithm for a NPGMNRs
initialize.model=function(x,y,k,method=NULL,true.init=NULL,p=1){
  n=length(y)
  BIC=1e6
  if(method=="1"){##Mixtures of Regression splines
    for(j in 1:1e3){
      m=list(BIC=1e6)
      try({m=mix.reg.splines(x,y,k)},silent=T)
      if(m$BIC<BIC){init.model=m$init.model0;BIC=m$BIC}
    }
  }
  if(method=="2"){##Mixtures of polynomial regressions
    for(j in 1:1e3){
      m=list(BIC=1e6)
      try({m=mix.poly(x,y,k,p)})
      if(m$BIC<BIC){init.model=m$init.model;BIC=m$BIC}
    }
  }
  if(method=="3"){##True values
    m0=true.init
    init.model=list(mu0=m0$mu,sigma20=m0$sigma2,pi0=m0$rho)
  }
  return(init.model)
}

##A function to fit a mixture of polynomial regressions (to initialize)
mix.poly=function(x,y,k,p){
  n=length(y)
  Xd=t(t(matrix(x,n,p+1))^(0:p))
  ##Initial state
  pi0=rep(1/k,k)
  sigma20=rgamma(k,1,1)^2
  Beta0=matrix(rnorm(k*(p+1)),p+1,k)
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))))))
  diff=1e6
  count=0
  while(diff>1e-10){
    ##E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-100)
    gn=g/rowSums(g)
    ##M-Step
    pi1=colSums(gn)/n
    Beta1=sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      Beta1=cbind(Beta1,solve(t(Xd)%*%W%*%Xd)%*%t(Xd)%*%W%*%y)
      sigma21=c(sigma21,sum(gn[,j]*(y-Xd%*%Beta1)^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-Xd%*%Beta1[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    Beta0=Beta1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==1e3) diff=1e-100
  }
  df=(p+1)*k+2*k-1
  BIC=-2*LogLik0+log(n)*df
  mu=Xd%*%Beta1
  r0=pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21));r=r0/rowSums(r0)
  if(which.max(pi1)==2) {pi1=pi1[2:1];sigma21=sigma21[2:1];Beta1=Beta1[,2:1]} 
  
  model0=list(r=r,Beta0=t(as.matrix(Beta1[-1,])),pi0=matrix(pi1,n,k,byrow=T),sigma20=matrix(sigma21,n,k,byrow=T),mu0=mu)
  return(list(init.model0=model0,BIC=BIC))
}

## A function to fit a Gaussian mixture model
GaussMix=function(y,x0,weights=NULL,model0){
  n=length(y)
  pi0=model0$pi0
  sigma20=model0$sigma20
  mu0=model0$mu0
  LL0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y,mu0[j],sqrt(sigma20)[j])))))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y,mu0[j],sqrt(sigma20)[j]))+1e-100;gn=g/rowSums(g)
    ##M-step
    pi1=colSums(weights*gn)/sum(weights)
    sigma21=mu1=NULL
    for(j in 1:k){
      W=weights*gn[,j]
      mu=sum(W*y)/sum(W)
      mu1=c(mu1,mu)
      sigma21=c(sigma21,sqrt(sum(W*(y-mu)^2)/sum(W)))
    }
    
    LL1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y,mu1[j],sqrt(sigma21)[j])))))
    diff=abs(LL1-LL0)
    LL0=LL1
    mu0=mu1
    pi0=pi1
    count=count+1
    sigma20=sigma21
    llk=c(llk,LL0)
    if(count==1e3) diff=1e-100
  }
  g=sapply(1:k,function(j) pi1[j]*dnorm(y,mu1[j],sqrt(sigma21)[j]));gn=g/rowSums(g)
  return(list(resp=gn,mix.mu=mu1,mix.sigma2=sigma21,mix.prop=pi1,LL=LL1))
}

##Roughness Approach
NPGMNRs_OB_EM=function(x,y,k,bw,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial State
  pi0=init.model0$pi0;
  mu0=init.model$mu0;
  sigma20=init.model$sigma20;
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  sigma20=sapply(1:k,function(j) approx(x,sigma20[,j],xgrid,rule=2)$y)
  pi0=sapply(1:k,function(j) approx(x,pi0[,j],xgrid,rule=2)$y)
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  R0=1e5
  llk=NULL
  while(diff>1e-10){
    ##E-Step
    zh=lapply(1:ngrid,function(t){
      model0=list(pi0=pi0[t,],sigma20=sigma20[t,],mu0=mu0[t,])
      model=GaussMix(y,xgrid[t],Kh[,t],model0);list(resp=model$resp,LL=model$LL)})
    ##M-Step
    mu=lapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        W=zh[[t]]$resp[,j]*Kh
        mh=colSums(W*y)/colSums(W)
      })
    })
    sigma2=lapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        W=zh[[t]]$resp[,j]*Kh
        mu2=approx(xgrid,mu[[t]][,j],xout=x,rule=2)$y
        sig2=colSums(W*(y-mu2)^2)/colSums(W)
      })
    })
    prop=lapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        W=zh[[t]]$resp[,j]*Kh
        prop=colSums(W)/colSums(Kh)
      })
    })
    ##Measuring the roughness
    for(t in 1:ngrid){
      R00=NULL
      for(j in 1:k) {R00=c(R00,Rough_curve(xgrid,mu[[t]][,j])$Rh)};R1=sum(log(R00))
      dist=max(apply(round(mu[[t]],4),1,function(x) abs(diff(x))))
      D=apply(mu[[t]],2,function(x) abs(trim(x)))
      id=apply(D,2,sum);id=which.min(id)
      ##Model space
      #No CRF can be a straight line (D>>0); two CRFs cannot be equal(dist>>0) 
      #How can I restrict the admission of straight lines? (sol:not all trimmed ranges should be zero)
      if(R1<R0 & dist>1e-1 & !all(D<0.5) & !all(D[,id]==0)){
        mu1=mu[[t]];pi1=prop[[t]];sigma21=sigma2[[t]];R0=R1;Rout=D
      }
    }
    ##Evaluating convergence
    mix.mu=mix.prop=mix.sigma2=NULL
    for(j in 1:k){
      mix.mu=cbind(mix.mu,approx(xgrid,mu1[,j],xout=x,rule=2)$y)
      mix.prop=cbind(mix.prop,approx(xgrid,pi1[,j],xout=x,rule=2)$y)
      mix.sigma2=cbind(mix.sigma2,approx(xgrid,sigma21[,j],xout=x,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2)[,j])))))
    diff=abs(LogLik1-LogLik0)
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    LogLik0=LogLik1
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  mu1=mix.mu
  ## Refitting to improve the two-stage estimates (mix.mu, mix.prop and mix.sigma2)
  out=refit(y,x,xgrid,k,mix.mu,mix.prop,mix.sigma2,bw);
  mu1=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  g=sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1)
  return(out)
}

##NaiveEM algorithm
NPGMNRs_Naive_EM=function(x,y,k,bw,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=init.model0$pi0
  mu0=init.model$mu0;
  sigma20=init.model$sigma20;
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  pi0=sapply(1:k,function(j) approx(x,pi0[,j],xgrid,rule=2)$y)
  sigma20=sapply(1:k,function(j) approx(x,sigma20[,j],xgrid,rule=2)$y)
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  tol=NULL
  while(diff>1e-10){
    ##local E-step
    g=lapply(1:ngrid,function(i){g=sapply(1:k,function(j) pi0[i,j]*dnorm(y-mu0[i,j],0,sqrt(sigma20[i,j])));gn=g/rowSums(g)})
    ##local M-step
    mu1=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]*Kh[t,]
        sum(w*y)/sum(w)
      })
    }))
    mu=sapply(1:k,function(j) approx(xgrid,mu1[,j],x,rule=2)$y)
    pi1=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]*Kh[t,]
        prop=sum(w)/sum(Kh[t,])
      })
    }))
    prop=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]*Kh[t,]
        sig2=sum(w*(y-mu)^2)/sum(w)
      })
    }))
    sigma2=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) prop[,j]*dnorm(y-mu[,j],0,sqrt(sigma2)[,j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e2) diff=1e-100
  }
  mu1=mu
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1)
  return(out)
}
           
NPGMNRs_Effective_EM=function(x,y,k,bw,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=init.model$pi0
  mu0=init.model$mu0;
  sigma20=init.model$sigma20
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])))))
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  count=0
  llk=NULL
  while(diff>1e-10){
    ##E-step
    g=sapply(1:k,function(j) pi0[,j]*dnorm(y-mu0[,j],0,sqrt(sigma20)[,j])+1e-100);gn=g/rowSums(g)
    ##M-step
    mu1=pi1=sigma21=NULL
    for(j in 1:k){
      W=gn[,j]*Kh
      pi1=cbind(pi1,colSums(W)/colSums(Kh))
      mh=colSums(W*y)/colSums(W)
      mu1=cbind(mu1,approx(xgrid,mh,x,rule=2)$y)
      sigma21=cbind(sigma21,colSums(W*(y-mu1[,j])^2)/colSums(W))
    }
    pi1=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=sapply(1:k,function(j) approx(xgrid,sigma21[,j],x,rule=2)$y)
    ##Evaluating convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[,j])))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    mu0=mu1
    pi0=pi1
    sigma20=sigma21
    count=count+1
    llk=c(llk,LogLik0)
    if(count==1e2) diff=1e-100
  }
  g=sapply(1:k,function(j) pi1[,j]*dnorm(y-mu1[,j],0,sqrt(sigma21[,j])));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1)
  return(out)
}
