##The following R package is required for fitting the Mixture-of-Experts (MOE) model
if(!require("mixtools")){install.packages("mixtools")}else{library(mixtools)}

###Kernel function
Kern<-function(x,x0,h){
  z=(x-x0)/h
  f=ifelse(abs(z)<=1, 0.75*(1 - z^2),0)/h
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

## A function to calculate the fitted BIC value
BIC=function(t,x,bw,K,LogLik){
  n=length(t)
  p=ncol(x)
  rk=2.1153;ck=0.045;
  dfm=(rk*abs(diff(range(t)))*ck)/bw
  dfn=(K-1)*dfm+(K+p*K)
  BIC=-2*LogLik+log(n)*dfn
  return(BIC)
}

##A one-step backfitting function
backfit=function(y,t,x,tgrid,k,Beta_init,pi_init,sigma2_init,bw,backfit=TRUE){
  n=length(y)
  z=as.matrix(x)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  Beta0=Beta_init
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  diff=1e10
  count=0
  while(diff>1e-10){
    #E-Step
    g=sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20[j])))
    gn=g/rowSums(g)
    #M-Step
    Beta1=sigma21=NULL
    for(j in 1:k){
      W=diag(gn[,j])
      Beta=solve(t(z)%*%W%*%z)%*%t(z)%*%W%*%y
      Beta1=cbind(Beta1,Beta)
      sigma21=c(sigma21,sum(gn[,j]*(y-z%*%Beta1[,j])^2)/sum(gn[,j]))
    }
    #Evaluate for convergence
    LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))))))
    diff=abs(LogLik0-LogLik1)
    LogLik0=LogLik1
    sigma20=sigma21
    Beta0=Beta1
    count=count+1
    if(count==5e2) diff=1e-100
  }
  #Re-estimating the non-parametric functions given the global parameter estimates
  mix.prop=pi0
  pi0=sapply(1:k,function(j) approx(t,pi0[,j],tgrid,rule=2)$y)
  if(backfit){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(tgrid,function(x0) Kern(t,x0,bw))
    ngrid=length(tgrid)
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      #gn=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi0[t,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      g=sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)
      ##local M-step
      pi1=sapply(1:k,function(j){
        #sapply(1:ngrid,function(t){
          W=gn[,j]*Kh
          #sum(W)/sum(Kh[,t])
          colSums(W)/colSums(Kh)
        })
      #})
      pi1=pi1/rowSums(pi1)
      prop=sapply(1:k,function(j) approx(tgrid,pi1[,j],t,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) prop[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      pi0=pi1
      count=count+1
      if(count==5e2) diff=1e-100
    }
  }
  out=list(Beta1=Beta1,pi1=mix.prop,sigma21=sigma21)
}

##An iterative backfitting function
backfit_fullyIter=function(y,t,x,tgrid,k,Beta_init,pi_init,sigma2_init,bw,backfit=TRUE){
  n=length(y)
  z=as.matrix(x)
  ngrid=length(tgrid)
  ##Estimating the global parameters given the non-parametric estimates
  pi0=pi_init
  sigma20=sigma2_init
  Beta0=Beta_init
  LogLikn=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  difff=1e6
  countN=0
  while(difff>1e-10){
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))  
    diff=1e10
    count=0
    while(diff>1e-10){
      #E-Step
      g=sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20[j])))
      gn=g/rowSums(g)
      #M-Step
      Beta1=sigma21=NULL
      for(j in 1:k){
        W=diag(gn[,j])
        Beta=solve(t(z)%*%W%*%z)%*%t(z)%*%W%*%y
        Beta1=cbind(Beta1,Beta)
        sigma21=c(sigma21,sum(gn[,j]*(y-z%*%Beta1[,j])^2)/sum(gn[,j]))
      }
      #Evaluate for convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      sigma20=sigma21
      Beta0=Beta1
      count=count+1
      if(count==5e2) diff=1e-100
    }
    
    #Re-estimating the non-parametric functions given the global parameter estimates
    mix.prop=pi0
    pi0=sapply(1:k,function(j) approx(t,pi0[,j],tgrid,rule=2)$y)
    LogLik0=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
    Kh=sapply(tgrid,function(x0) Kern(t,x0,bw))
    diff=1e10
    count=0
    while(diff>1e-10){
      ##local E-step
      g=lapply(1:ngrid,function(t){g=sapply(1:k,function(j) pi0[t,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);gn=g/rowSums(g)})
      #g=sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21[j]))+1e-100);
      #gn=g/rowSums(g)
      ##local M-step
      pi1=sapply(1:k,function(j){
        sapply(1:ngrid,function(t){
          W=g[[t]][,j]*Kh[,t]
          sum(W)/sum(Kh[,t])
          #colSums(W)/colSums(Kh[,t])
        })
      })
      pi1=pi1/rowSums(pi1)
      prop=sapply(1:k,function(j) approx(tgrid,pi1[,j],t,rule=2)$y)
      ##Evaluating convergence
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) prop[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j])))))
      diff=abs(LogLik0-LogLik1)
      LogLik0=LogLik1
      pi0=pi1
      count=count+1
      if(count==5e2) diff=1e-100
    }
    pi1=sapply(1:k,function(j) approx(tgrid,pi1[,j],t,rule=2)$y)
    ##Evaluate for overall convergence
    difff=abs(LogLik1-LogLikn)
    LogLikn=LogLik1
    Beta0=Beta1
    pi0=pi1
    sigma20=sigma21
    countN=countN+1
    if(countN==1e2) difff=1e-100
  }
  out=list(Beta1=Beta1,pi1=mix.prop,sigma21=sigma21)
}

##A function to initialize the fitting algorithm in the function SPGMRVPs_MB_ECM
initialize.model=function(x,t,y,k,method=NULL,true.init=NULL,p=1){
  n=length(y)
  BIC=1e6
  if(method=="1"){##Mixtures-of-Experts
    for(j in 1:1e2){
      m=list(BIC=1e6)
      try({m=HME(y,x,t,k=2)})
      init.model=list(mu0=m$mix.mu,sigma20=m$mix.sigma2,pi0=m$mix.prop,Beta0=m$mix.beta)
      if(m$BIC<BIC){init.model=init.model;BIC=m$BIC}
    }
  }
  if(method=="2"){##Mixtures of polynomial regressions
    for(j in 1:1e2){
      m=list(BIC=1e6)
      try({m=mix.poly(x,y,k,p)})
      if(m$BIC<BIC){init.model=m$init.model;BIC=m$BIC}
    }
    init.model$pi0=matrix(init.model$pi0,n,k,byrow=T)
  }
  if(method=="2"){##True values
    m0=true.init
    init.model=list(mu0=m0$mix.mu,sigma20=m0$mix.sigma2,pi0=m0$mix.prop,Beta0=m0$mix.beta)
  }
  return(init.model)
}

##A function to fit the mixture of polynomial regressions (to initialize)
mix.poly=function(x,y,k,p){
  n=length(y)
  xd=poly(x,p)
  fit=regmixEM(y,x,k=k)
  pi1=fit$lambda;sigma21=fit$sigma^2;Beta1=fit$beta
  df=(p+1)*k+2*k-1
  BIC=-2*fit$loglik+log(n)*df
  mu=(fit$x)%*%Beta1
  r0=sapply(1:k,function(j) pi1[j]*dnorm(y-mu[,j],0,sqrt(sigma21)[j]));r=r0/rowSums(r0)
  if(which.max(pi1)==2) {pi1=pi1[2:1];sigma21=sigma21[2:1];Beta1=Beta1[,2:1]} 
  
  model0=list(r=r,Beta0=Beta1,pi0=pi1,sigma20=sigma21,mu0=mu)
  return(list(init.model0=model0,BIC=BIC))
}

## A modified function to fit the Mixture-of-Experts using the hmeEM function from 'mixtools' package
HME=function(y,x,t,k){
  n=length(y)
  model=hmeEM(y,x,k=k);
  w=as.matrix(model$w)
  z=cbind(rep(1,n),x)
  mix.prop=1/(1+exp(-z%*%w));mix.prop=cbind(mix.prop,1-mix.prop)
  #mix.prop=sapply(2:k,function(j) 1/(1+exp(-z%*%w[,j])));mix.prop=cbind(mix.prop,1-mix.prop)
  mix.beta=model$beta#model.param[1:2,]
  mix.sigma2=(model$sigma)^2#model.param[3,]^2
  res=order(mix.sigma2);mix.beta=mix.beta[,res];mix.prop=mix.prop[,res];mix.sigma2=mix.sigma2[res]
  mu=sapply(1:k,function(j) z%*%mix.beta[,j])
  g=sapply(1:k,function(j) mix.prop[,j]*dnorm(y-z%*%mix.beta[,j],0,sqrt(mix.sigma2)[j]));gn=g/rowSums(g)
  LL=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-z%*%mix.beta[,j],0,sqrt(mix.sigma2)[j])))))
  df=(k+ncol(z)*k)+(k-1)*nrow(w)
  BIC=-2*LL+log(n)*df
  return(list(resp=gn,mix.prop=mix.prop,mix.mu=mu,mix.beta=mix.beta,mix.sigma2=mix.sigma2,LL=LL,BIC=BIC))
}

## A function to compute the conditional distribution for a GMLRs model
GMLR=function(y,z,mix.beta,mix.prop,mix.sigma){
  k=length(mix.prop)
  out=rowSums(sapply(1:k,function(j) mix.prop[j]*dnorm(y-z%*%mix.beta[,j],0,mix.sigma[j])))
  return(out)
}

##A function to fit the SPGMRVPs using the model-based approach via the ECM algorithm
SPGMRVPs_MB_ECM=function(x,t,y,k,bw,tgrid,init.model,lmd_0=1e-5){
  n=length(y)
  x=as.matrix(x)
  z=cbind(rep(1,n),x)
  p=ncol(z)
  u=t
  ngrid=length(tgrid)
  grid=1:ngrid
  tgrid0=tgrid
  ##Initial State
  lmd0=rep(1/ngrid,ngrid)
  pi0=init.model$pi0
  Beta0=init.model$Beta0;
  sigma20=init.model$sigma20
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  mix.Beta0=lapply(1:k,function(j) matrix(Beta0[,j],n,p,byrow=T))
  Beta0=lapply(1:k,function(j) matrix(Beta0[,j],ngrid,p,byrow=T))
  pi0=sapply(1:k, function(j) approx(u,pi0[,j],tgrid,rule=2)$y)
  sigma20=matrix(sigma20,ngrid,k,byrow = T)
  Kh=sapply(tgrid,function(u0) Kern(u,u0,bw))
  Kh0=Kh
  diff=1e6
  count=0
  llk=NULL
  tol=NULL
  while(diff>1e-10){
    ##E-Step
    v=sapply(1:length(grid),function(t){
      beta0=sapply(1:k,function(j) Beta0[[j]][t,])
      lmd0[t]*GMLR(y,z,beta0,pi0[t,],sqrt(sigma20)[t,])+1e-300});vh=v/rowSums(v)
    zh=lapply(1:length(grid),function(t){
      beta0=sapply(1:k,function(j) Beta0[[j]][t,])
      g=sapply(1:k,function(j) pi0[t,j]*dnorm(y-z%*%beta0[,j],0,sqrt(sigma20[t,j])))+1e-300
      ;gh=g/rowSums(g)})
                   
    ##M-Step
    lmd1=colSums(vh)/n
    grid=which(lmd1>=lmd_0)
    lmd1=lmd1[grid]
    tgrid=tgrid[grid]
    Kh=Kh[,grid]
                   
    ### Beta
    Beta=lapply(1:k,function(j){
        #Next-level (a list with size=k within each top-level list)
        b=t(sapply(1:length(grid),function(t){
          id=grid[t]
          W=diag(lmd1[t]*zh[[id]][,j]*Kh[,t])
          solve(t(z)%*%W%*%z)%*%t(z)%*%W%*%y
        }))
      })
    ### sigma2
    sigma2=sapply(1:k,function(j){
    if(p>1){
      eta=apply((Beta)[[j]],2,function(x) approx(tgrid,x,xout=u,rule=2)$y)
    }else{
      eta=approx(tgrid,(Beta)[[j]],xout=u,rule=2)$y
    }
    sapply(1:length(grid),function(t){
    id=grid[t]
    W=lmd1[t]*zh[[id]][,j]*Kh[,t]
     yh=rowSums(z*eta)
     res2=c((y-yh)^2)
    sig2=sum(W*res2)/sum(W)
    })
  })

  ### mix proportion
  prop=sapply(1:k,function(j){
    sapply(1:length(grid),function(t){
      id=grid[t]
      W=lmd1[t]*zh[[id]][,j]*Kh[,t]
      rho=sum(W)/sum(lmd1[t]*Kh[,t])
    })
  })
  prop=prop/rowSums(prop)
  ##Evaluating convergence
  if(p>1){
    mix.beta=lapply(1:k,function(j) apply(Beta[[j]],2,function(x) approx(tgrid,x,xout=u,rule=2)$y))
  }else{
    mix.beta=lapply(1:k,function(j) approx(tgrid,Beta[[j]],xout=u,rule=2)$y)
  }
  mix.prop=mix.sigma2=NULL
  for(j in 1:k){
    mix.prop=cbind(mix.prop,approx(tgrid,prop[,j],xout=u,rule=2)$y)
    mix.sigma2=cbind(mix.sigma2,approx(tgrid,sigma2[,j],xout=u,rule=2)$y)
  }
  LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-rowSums(z*mix.beta[[j]]),0,sqrt(mix.sigma2)[,j])))))
  diff=abs(LogLik1-LogLik0)
  dif=max(abs(lmd0[grid]-lmd1))
  lmd0=lmd1
  sigma20=sigma2
  pi0=prop
  Beta0=Beta
  LogLik0=LogLik1
  count=count+1
  tol=c(tol,LogLik1)
  if(count==1e2) diff=1e-100
  }
  beta1=sapply(1:k,function(j) colMeans(mix.beta[[j]]))
  out=backfit(y,u,z,tgrid0,k,beta1,mix.prop,colMeans(mix.sigma2),bw);
  Beta1=out$Beta1;pi1=out$pi1;sigma21=out$sigma21;
  res=order(sigma21);Beta1=Beta1[,res];pi1=pi1[,res];sigma21=sigma21[res]
  g=sapply(1:k,function(j) pi1[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  BIC=BIC(t,cbind(1,x),bw,k,LL1)
  mu=sapply(1:k,function(j) z%*%Beta1[,j])
  out=list(resp=gn,mix.mu=mu,mix.prop=pi1,mix.beta=Beta1,mix.sigma2=sigma21,LL=LL1,BIC=BIC)
  return(out)
}

##A function to fit the SPGMRVPs using the Naive EM algorithm
SPGMRVPs_Naive_EM=function(x,t,y,k,bw,tgrid,init.model){
  n=length(y)
  x=as.matrix(x)
  z=cbind(rep(1,n),x)
  p=ncol(z)
  u=t
  ngrid=length(tgrid)
  ##Initial State
  pi0=init.model$pi0
  Beta0=init.model$Beta0;
  sigma20=init.model$sigma20
  ##
  LogLik0=sum(log(rowSums(sapply(1:k,function(j) pi0[,j]*dnorm(y-z%*%Beta0[,j],0,sqrt(sigma20)[j])))))
  mix.Beta0=lapply(1:k,function(j) matrix(Beta0[,j],n,p,byrow=T))
  Beta0=lapply(1:k,function(j) matrix(Beta0[,j],ngrid,p,byrow=T))
  pi0=sapply(1:k, function(j) approx(u,pi0[,j],tgrid,rule=2)$y)
  sigma20=matrix(sigma20,ngrid,k,byrow = T)
  Kh=sapply(tgrid,function(u0) Kern(u,u0,bw))
  Kh0=Kh
  diff=1e6
  count=0
  llk=NULL
  tol=NULL
  while(diff>1e-10){
    ##E-Step
      gn=lapply(1:ngrid,function(t){
        beta0=sapply(1:k,function(j) Beta0[[j]][t,])
        g=sapply(1:k,function(j) pi0[t,j]*dnorm(y-z%*%beta0[,j],0,sqrt(sigma20[t,j])))+1e-300
        ;gh=g/rowSums(g)})
      ##M-Step
      ##Beta
      Beta=lapply(1:k,function(j){
        #Next-level (a list with size=k within each top-level list)
        b=t(sapply(1:ngrid,function(t){
          W=diag(gn[[t]][,j]*Kh[,t])
          solve(t(z)%*%W%*%z)%*%t(z)%*%W%*%y
        }))
      })
      ##sigma2
      sigma2=sapply(1:k,function(j){
        if(p>1){
          eta=apply((Beta)[[j]],2,function(x) approx(tgrid,x,xout=u,rule=2)$y)
        }else{
          eta=approx(tgrid,(Beta)[[j]],xout=u,rule=2)$y
        }
        sapply(1:ngrid,function(t){
          W=gn[[t]][,j]*Kh[,t]
          yh=rowSums(z*eta)
          res2=c((y-yh)^2)
          sig2=sum(W*res2)/sum(W)
        })
      })
      ##mixing proportion
      prop=sapply(1:k,function(j){
        sapply(1:ngrid,function(t){
          W=gn[[t]][,j]*Kh[,t]
          rho=sum(W)/sum(Kh[,t])
        })
      })
      prop=prop/rowSums(prop)
      ##Evaluating convergence
      if(p>1){
        mix.beta=lapply(1:k,function(j) apply(Beta[[j]],2,function(x) approx(tgrid,x,xout=u,rule=2)$y))
      }else{
        mix.beta=lapply(1:k,function(j) approx(tgrid,Beta[[j]],xout=u,rule=2)$y)
      }
      mix.prop=mix.sigma2=NULL
      for(j in 1:k){
        mix.prop=cbind(mix.prop,approx(tgrid,prop[,j],xout=u,rule=2)$y)
        mix.sigma2=cbind(mix.sigma2,approx(tgrid,sigma2[,j],xout=u,rule=2)$y)
      }
      LogLik1=sum(log(rowSums(sapply(1:k,function(j) mix.prop[,j]*dnorm(y-rowSums(z*mix.beta[[j]]),0,sqrt(mix.sigma2)[,j])))))
      diff=abs(LogLik1-LogLik0)
      sigma20=sigma2
      pi0=prop
      Beta0=Beta
      LogLik0=LogLik1
      count=count+1
      tol=c(tol,LogLik1)
      if(count==1e2) diff=1e-100
  }
  beta1=sapply(1:k,function(j) colMeans(mix.beta[[j]]))
  out=backfit(y,u,z,tgrid,k,beta1,mix.prop,colMeans(mix.sigma2),bw,backfit=F);
  Beta1=out$Beta1;pi1=out$pi1;sigma21=out$sigma21;
  res=Dist(pi1);Beta1=Beta1[,res$id];pi1=pi1[,res$id];sigma21=sigma21[res$id]
  g=sapply(1:k,function(j) pi1[,j]*dnorm(y-z%*%Beta1[,j],0,sqrt(sigma21)[j]));gn=g/rowSums(g)
  LL1=sum(log(rowSums(g)))
  BIC=BIC(t,cbind(1,x),bw,k,LL1)
  mu=sapply(1:k,function(j) z%*%Beta1[,j])
  out=list(resp=gn,mix.mu=mu,mix.prop=pi1,mix.beta=Beta1,mix.sigma2=sigma21,LL=LL1,BIC=BIC)
  return(out)
}
