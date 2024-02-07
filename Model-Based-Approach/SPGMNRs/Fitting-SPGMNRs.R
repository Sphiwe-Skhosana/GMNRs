if(!require("mixtools")){install.packages("mixtools")}else{library(mixtools)}
if(!require("locpol")){install.packages("locpol")}else{library(locpol)}
library(splines)
library(MASS)
###Kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  f=dnorm(z)
  out=f
  if(sum(f)>0){out=f/sum(f)};
  return(out)
}

##Normalizing functions

Normalize=function(x){
  x/sum(x)
}

minmaxscaler=function(x){
  (x-min(x))/(max(x)-min(x))
}

standardizer=function(x){
  z=(x-mean(x))/sd(x)
  return(z)
}

##Conditional distribution of y|x
cond_dist=function(y,mix.mu,mix.prop,mix.sigma2){
  k=length(mix.prop)
  rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2[j]))))
}

##A function to calculate the BIC
BIC=function(x,bw,K,LogLik,dfn){
  n=length(x)
  rk=2.5375;ck=0.7737;
  dfm=(rk*abs(diff(range(x)))*ck)/bw
  dfp=2*K-1
  df=K*dfm+2*K-1
  BIC1=-2*LogLik+log(n)*df
  BIC2=-2*LogLik+log(n)*dfn
  return(c(BIC1,BIC2))
}

###Kernel polynomial smoothing function
polynomial.smoother=function(x,y,xgrid,d,W){
  ##d: denotes the degree of the polynomial
  n=length(x);ngrid=length(xgrid)
  fit=sapply(1:ngrid,function(i){
    x0=xgrid[i];W=diag(W[,i])+1e-100*diag(n)
    z=(x-x0);Z=t(t(matrix(z,n,d+1))^(0:d))
    (solve(t(Z)%*%W%*%Z)%*%t(Z)%*%W%*%y)[1]
  })
  return(fit)
}

###Kernel polynomial smoother matrix
polynomial.smoother.matrix=function(x,x.grid,d,W){
  n=length(x)
  S=sapply(1:length(x.grid),function(i){
    w=diag(W[,i]);
    X=t(t(matrix((x-x.grid[i]),n,d+1))^(0:d))
    (solve(t(X)%*%w%*%X)%*%t(X)%*%w)[1,]
  })
  return(S)
}

##Backfitting function
backfit=function(y,x,xgrid,d,k,mh,pi_init,sigma2_init,bw,backfit=TRUE){
  n=length(y)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20)[j])))))
  diff=1e10
  count=0
  while(diff>1e-10){
    #E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20[j]))+1e-300)
    gn=g/rowSums(g)
    #M-Step
    pi1=colSums(gn)/n
    sigma21=NULL
    for(j in 1:k){
      sigma21=c(sigma21,sum(gn[,j]*(y-mh[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==5e2) diff=1e-100
  }
  
  #Re-estimating the non-parametric functions given the global parameter estimates
  ##Initialiaze the algorithm using the responsibilities from the previous stage
  ##Since the parameter estimates are well-labelled
  mu=mh
  mu0=mh
  #mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  if(backfit){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
    ngrid=length(xgrid)
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      #gn=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[t,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      #mu1=t(sapply(1:ngrid,function(t){
      mu1=sapply(1:k,function(j){
        W=gn[,j]*Kh
        mh=colSums(W*y)/colSums(W)
        #mh=sum(W*y)/sum(W)
        #approx(xgrid,mh,xout=x,rule=2)$y
      })
      #}))
      mu=sapply(1:k,function(j) approx(xgrid,mu1[,j],x,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mu
      count=count+1
      if(count==5e2) diff=1e-100
    }
  }
  out=list(mu=mu,pi1=pi1,sigma21=sigma21)
}

##An iterative backfitting function
backfit_fullyIter=function(y,x,xgrid,d,k,mh,pi_init,sigma2_init,bw){
  n=length(y)
  ngrid=length(xgrid)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  LogLikn=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20)[j])))))
  difff=1e6
  countN=0
  while(difff>1e-10){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20)[j])))))  
    diff=1e10
    count=0
    while(diff>1e-10){
      #E-Step
      g=sapply(1:k,function(j) pi0[j]*dnorm(y-mh[,j],0,sqrt(sigma20[j]))+1e-300)
      gn=g/rowSums(g)
      #M-Step
      pi1=colSums(gn)/n
      sigma21=NULL
      for(j in 1:k){
        sigma21=c(sigma21,sum(gn[,j]*(y-mh[,j])^2)/sum(gn[,j]))
      }
      #Evaluate for convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mh[,j],0,sqrt(sigma21[j]))))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      sigma20=sigma21
      pi0=pi1
      count=count+1
      if(count==5e2) diff=1e-100
    }
    
    #Re-estimating the non-parametric functions given the global parameter estimates
    mu0=mh
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j],0,sqrt(sigma21)[j])))))
    mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
    Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
    diff=1e10
    count=0
    while(diff>1e-10){
      #local E-step
      gn=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[t,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      #g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu0[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      mugrid=t(sapply(1:ngrid,function(t){
        mugrid=sapply(1:k,function(j){
          W=gn[[t]][,j]*Kh
          #mh=colSums(W*y)/colSums(W)
          mh=sum(W*y)/sum(W)
          #approx(xgrid,mh,xout=x,rule=2)$y
        })
        }))
      mu=sapply(1:k,function(j) approx(xgrid,mugrid[,j],x,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      mu0=mugrid
      count=count+1
      if(count==5e2) diff=1e-100
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

##A function to perform local polynomial smoothing or estimation
local.polynomial.smoother=function(x,y,xgrid,bw,d,W){
  library(locpol)
  n=length(y)
  g=locPolSmootherC(x,y,xgrid,bw,d,gaussK,weig=W)
  return(g)
}

###A function to compute the roughness of a function
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

###Initialization function
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
    for(j in 1:1e2){
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

##A function to fit the mixture of polynomial regressions (to initialize)
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
  
  model0=list(r=r,Beta0=Beta1,pi0=pi1,sigma20=sigma21,mu0=mu)
  return(list(init.model0=model0,BIC=BIC))
}

###A function to fit the mixture of regression splines (to initialize)
mix.reg.splines=function(x,y,k){
  n=length(y)
  Xd=bs(x,knots=quantile(unique(x),probs=seq(0.2,0.8,0.2)),intercept=T)##B-spline basis matrix
  d=ncol(Xd);
  ##Initial state
  pi0=rep(1/k,k)
  sigma20=rgamma(k,1,1)^2
  Beta0=matrix(rnorm(k*d),d,k)
  LogLik0=sum(log(rowSums(pi0*dnorm(y-Xd%*%Beta0,0,sqrt(sigma20)))))
  diff=1e6
  ll=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    g=sapply(1:k,function(j) pi0[j]*dnorm(y-Xd%*%Beta0[,j],0,sqrt(sigma20[j]))+1e-300)
    gn=g/rowSums(g)
    ##M-Step
    pi1=colSums(gn)/n
    Beta1=sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      Beta1=cbind(Beta1,solve(t(Xd)%*%W%*%Xd)%*%t(Xd)%*%W%*%y)
      sigma21=c(sigma21,sum(gn[,j]*(y-Xd%*%Beta1[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21)))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    Beta0=Beta1
    sigma20=sigma21
    pi0=pi1
    count=count+1
    if(count==1e2) diff=1e-100
  }
  df=k*d+2*k-1
  df_reg=k*d
  BIC=-2*LogLik0+df*log(n)
  mu1=Xd%*%Beta1
  r0=pi1*dnorm(y-Xd%*%Beta1,0,sqrt(sigma21));r=r0/rowSums(r0)
  R2=Rsquared(y,mu1,r)
  if(which.max(pi1)==2) {pi1=pi1[2:1];sigma21=sigma21[2:1];Beta1=Beta1[,2:1]} 
  model0=list(r=r,Beta0=matrix(Beta1,d,k),pi0=pi1,sigma20=sigma21,mu0=mu1,LL=LogLik1,R2=R2,BIC=BIC)
  return(list(r=r,BIC=BIC,mu=mu1,Beta=Beta1,pi=pi1,sigma2=sigma21,init.model0=model0,df_reg=df_reg))
}

GMM=function(y,mix.mu,mix.prop,mix.sigma){
  k=length(mix.mu)
  out=rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-mix.mu[j],0,mix.sigma[j])))
  return(out)
}

###GMLRs
GaussLinMix=function(x,y,k,init.model0){
  n=length(y)
  prop0=init.model0$pi0;sigma0=sqrt(init.model0$sigma20);beta0=init.model0$Beta0[1:2,]
  fit=regmixEM(y,x,lambda=prop0,sigma=sigma0,beta=beta0)
  mu1=(fit$x)%*%(fit$beta);beta=(fit$beta);pi1=fit$lambda;sigma21=fit$sigma^2
  res=Dist(mu1);mu1=mu1[,res$id];beta1=beta[,res$id];pi1=pi1[res$id];sigma21=sigma21[res$id]
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  R2=Rsquared(y,mu1,gn)$RMSE
  df=4*k-1
  LL=fit$loglik
  BIC=-2*LL+df*log(n)
  out=list(resp=gn,mix.prop=pi1,mix.sigma2=sigma21,mix.mu=mu1,beta=beta1,x=fit$x,LL=fit$loglik,R2=R2,df=df,BIC=BIC)
  return(out)
}

###A function to fit the SPGMNRs model using the MB-EM algorithm
SPGMNR_MB_EM=function(x,y,k,bw,d,xgrid,init.model,lmd_0=NULL){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  lmd0=rep(1/ngrid,ngrid)
  pi0=init.model$pi0;pi0=matrix(pi0,ngrid,k,byrow=T)
  mu0=init.model$mu0;
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  sigma20=init.model$sigma20;sigma20=matrix(sigma20,ngrid,k,byrow=T)
  ##
  LogLik0=sum(log(rowSums(sapply(1:ngrid,function(t) lmd0[t]*rowSums(sapply(1:k,function(j) pi0[t,j]*dnorm(y-mu0[t,j],0,sqrt(sigma20)[t,j])))))))
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=dif=1e6
  count=0
  tol=NULL
  while(diff>1e-10){
    ##E-Step
    v=sapply(1:ngrid,function(t) lmd0[t]*GMM(y,mu0[t,],pi0[t,],sqrt(sigma20)[t,]))+1e-300;vh=v/rowSums(v)
    zh=lapply(1:ngrid,function(t){
      g=sapply(1:k,function(j) pi0[t,j]*dnorm(y,mu0[t,j],sqrt(sigma20[t,j])))+1e-300
      ;gh=g/rowSums(g)})
    ##M-Step
    lmd1=colSums(vh)/n
    ##mu
    mu=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]
        mh=local.polynomial.smoother(x,y,xgrid[t],bw,d,W)[,2]
        #sum(W*y)/sum(W)
      })
    })
    ##sigma2
    sigma2=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]*Kh[,t]
        yh=mu[t,j]
        res2=c((y-yh)^2)
        sig2=sum(W*res2)/sum(W)
      })
    })
    
    prop=sapply(1:k,function(j){
      sapply(1:ngrid,function(t){
        W=lmd1[t]*zh[[t]][,j]*Kh[,t]
        rho=sum(W)/sum(lmd1[t]*Kh[,t])
      })
    })
    prop=prop/rowSums(prop)
    ##Evaluating convergence
    mix.mu=mix.prop=mix.sigma2=NULL
    for(j in 1:k){
      mix.prop=cbind(mix.prop,approx(xgrid,prop[,j],xout=x,rule=2)$y)
      mix.sigma2=cbind(mix.sigma2,approx(xgrid,sigma2[,j],xout=x,rule=2)$y)
      mix.mu=cbind(mix.mu,approx(xgrid,mu[,j],xout=x,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2)[,j])))))
    diff=abs(LogLik1-LogLik0)
    dif=max(abs(lmd0-lmd1))
    lmd0=lmd1
    sigma20=sigma2
    pi0=prop
    mu0=mu
    LogLik0=LogLik1
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e2) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi0),colMeans(sigma20),bw);
  mu2=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.mu1=mu1,mix.sigma2=sigma21,df=df,BIC=BIC[2],LL=LL1)
  return(out)
}

###A function to fit SPGMNRs using the MB-ECM algorithm
SPGMNRs_MB_ECM=function(x,y,k,bw,d,xgrid,init.model,lmd_0=1e-5){
  n=length(y)
  ngrid=length(xgrid)
  grid=1:ngrid
  xgrid0=xgrid
  ##Initial state
  lmd0=rep(1/ngrid,ngrid)
  pi0=init.model$pi0;pi0=matrix(pi0,ngrid,k,byrow=T)
  mu0=init.model$mu0;
  mu0=sapply(1:k,function(j) approx(x,mu0[,j],xgrid,rule=2)$y)
  sigma20=init.model$sigma20;sigma20=matrix(sigma20,ngrid,k,byrow=T)
  ##
  LogLik0=sum(log(rowSums(sapply(1:ngrid,function(t) lmd0[t]*rowSums(sapply(1:k,function(j) pi0[t,j]*dnorm(y-mu0[t,j],0,sqrt(sigma20)[t,j])))))))
  Kh=sapply(xgrid,function(x0) Kern(x,x0,bw))
  diff=1e6
  tol=NULL
  count=0
  while(diff>1e-10){
    ##E-Step
    v=sapply(1:length(grid),function(t) lmd0[t]*GMM(y,mu0[t,],pi0[t,],sqrt(sigma20)[t,]))+1e-300;vh=v/rowSums(v)
    zh=lapply(1:length(grid),function(t){
      g=sapply(1:k,function(j) pi0[t,j]*dnorm(y,mu0[t,j],sqrt(sigma20[t,j])))+1e-300
      ;gh=g/rowSums(g)})
    ##M-Step
    lmd1=colSums(vh)/n
    grid=which(lmd1>=lmd_0)
    lmd1=lmd1[grid]
    xgrid=xgrid[grid]
    Kh=Kh[,grid]
    ##mu
    mu=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=lmd1[t]*zh[[id]][,j]*Kh[,t]
        mh=local.polynomial.smoother(x,y,xgrid[t],bw,d,W)[,2]
        #sum(W*y)/sum(W)
      })
    })
    ##sigma2
    sigma2=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=lmd1[t]*zh[[id]][,j]*Kh[,t]
        yh=mu[t,j]
        res2=c((y-yh)^2)
        sig2=sum(W*res2)/sum(W)
      })
    })
    
    prop=sapply(1:k,function(j){
      sapply(1:length(grid),function(t){
        id=grid[t]
        W=zh[[id]][,j]*Kh[,t]
        rho=sum(W)/sum(lmd1[t]*Kh[,t])
      })
    })
    prop=prop/rowSums(prop)
    ##Evaluating convergence
    mix.mu=mix.prop=mix.sigma2=NULL
    for(j in 1:k){
      mix.prop=cbind(mix.prop,approx(xgrid,prop[,j],xout=x,rule=2)$y)
      mix.sigma2=cbind(mix.sigma2,approx(xgrid,sigma2[,j],xout=x,rule=2)$y)
      mix.mu=cbind(mix.mu,approx(xgrid,mu[,j],xout=x,rule=2)$y)
    }
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-mix.mu[,j],0,sqrt(mix.sigma2)[,j])))))
    diff=abs(LogLik1-LogLik0)
    dif=max(abs(lmd0[grid]-lmd1))
    lmd0=lmd1
    sigma20=sigma2
    pi0=prop
    mu0=mu
    LogLik0=LogLik1
    count=count+1
    tol=c(tol,LogLik1)
    if(count==1e3) diff=1e-100
  }
  mu1=mix.mu
  out=backfit(y,x,xgrid0,d,k,mu1,colMeans(pi0),colMeans(sigma20),bw);
  mu2=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu2[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu2,mix.sigma2=sigma21,df=df,BIC=BIC[2],LL=LL1)
  return(out)
}

##A function to fit the SPGMNRs model using the local EM (LEM) algorithm 
SPGMNRs_LEM=function(x,y,k,bw,d,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=t((init.model$pi0)*t(matrix(1,n,k)))
  mu0=init.model$mu0;
  sigma20=t((init.model$sigma20)*t(matrix(1,n,k)))
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
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi1),colMeans(sigma21),bw);
  mu1=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21[j])));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  Kh=sapply(x,function(x0) Kern(x,x0,bw))
  df=sum(sapply(1:k,function(j){w=gn[,j]*Kh;S=polynomial.smoother.matrix(x,x,d,w);sum(diag(S))}))
  BIC=BIC(x,bw,k,LL1,df)
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,df=df,LL=LL1,BIC=BIC)
  return(out)
}

##A function to fit the SPGMNRs using the Naive EM algorithm
SPGMNRs_Naive_EM=function(x,y,k,bw,d,xgrid,init.model){
  n=length(y)
  ngrid=length(xgrid)
  ##Initial state
  pi0=init.model0$pi0;pi0=matrix(pi0,n,k,byrow=T)
  mu0=init.model$mu0;
  sigma20=init.model$sigma20;sigma20=matrix(sigma20,n,k,byrow=T)
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
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,y,xgrid[t],bw,0,w)[,2]
      })
    }))
    mu=sapply(1:k,function(j) approx(xgrid,mu1[,j],x,rule=2)$y)
    pi1=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,rep(1,n),xgrid[t],bw,0,w)[,2]
      })
    }))
    prop=sapply(1:k,function(j) approx(xgrid,pi1[,j],x,rule=2)$y)
    sigma21=t(sapply(1:ngrid,function(t){
      sapply(1:k,function(j){
        w=g[[t]][,j]
        mh=local.polynomial.smoother(x,(y-mu1[t,j])^2,xgrid[t],bw,0,w)[,2]
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
  out=backfit(y,x,xgrid,d,k,mu1,colMeans(pi1),colMeans(sigma21),bw,backfit = F);
  mu1=out$mu;pi1=out$pi1;sigma21=out$sigma21;
  g=sapply(1:k,function(j) pi1[j]*dnorm(y-mu1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  out=list(resp=gn,mix.prop=pi1,mix.mu=mu1,mix.sigma2=sigma21,LL=LL1)
  return(out)
}
