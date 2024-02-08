# Introduction
This folder contains all the ```R``` code used to fit semi-parametric mixtures of regressions with varying mixing proportions (SPGMRVPs) uisng the proposed model-based approach, the Mixtures-of-Experts via the parametric EM algorithm (Benanglia et. al., 2010) and the naive EM algorithm.

# Description of the code
In this section, we provide a brief description of the code (contained in the ```R``` script Fitting-SPGMNRs.R):
* The ```SPGMRVPs_MB_ECM(x,t,y,k,bw,xgrid,init.model,lmd=1e-5)``` function fits the SPGMNRs model using the model-based ECM approach.
* The ```SPGMRVPs_Naive_EM (x,t,y,k,bw,xgrid,init.model)``` function fits the SPGMNRs model using the naive EM algorithm without consideration to the label-switching problem.
* The ```HME (y,x,t,k)``` function fits the Mixture-of-Experts model via the parametric EM algorithm using the ```hmeEM``` function in the ```R``` package ```mixtools``` (. The LEM algorithm makes use of the effective EM algorithm to account for the label-switching problem.

  #### Arguments (inputs)
  + ```x``` is a vector of length $n$ that consists of the covariate values. Note that the function can only take only one covariate.
  + ```y``` is a vector of length $n$ that consists of the response variable values
  + ```k``` is the number of components
  + ```bw``` is the bandwidth 
  + ```xgrid``` is a vector of length $N\leq n$ that consists of a set of local grid points
  + ```init.model``` is a list object that contains the model to initialize the algorithm. It is obtained as the output of the ```initialize.model``` function. See the description of the ```initialize.model``` function below.
  + ```lmd``` specifies the threshold parameter value. The default is ```lmd=1e-5```.
    
  #### Values (outputs)
  All of the above functions return a list with the following items:
  + ```resp``` is an $n\times k$ matrix of the fitted responsibilities
  + ```mix.prop``` is the fitted mixing proportions
  + ```mix.mu``` is the fitted component regression functions (CRFs)
  + ```mix.sigma2``` is the fitted variances
  + ```df``` is the total degrees of freedom for the fitted CRFs
  + ```BIC``` is the value of the fitted Bayesian information criterion
  + ```LL``` is the fitted log-likelihood
* The ```initialize.model(x,y,k,method=1,true.init.model=NULL,p=1)``` function fits/computes a mixture of regressions model using one of three methods specified by the argument ```method```. The fitted model is used to initialize any of the above functions used to fit the SPGMNRs model.
    #### Arguments (inputs)
  + ```x``` is a vector of length $n$ that consists of the covariate values. Note that the function can only take only one covariate.
  + ```y``` is a vector of length $n$ that consists of the response variable values
  + ```k``` is an integer that specifies the number of components
  + ```method``` is an integer that specifies the initialization method that should be used. The choices are ```1- mixture of regression splines```, ```2 - mixture of polynomial regressions``` and ```3 - the true SPGMNRs model (if known)```
  + ```true.init.model``` is a list object that contains the true model
  + ```p``` If ```method=2``` then this parameter specifies the degree of the polynomial.    
  #### Values (output)
  The ```initialize.model``` function returns a list object ```init.model0``` with the initial fitted model.
# References
1. Xiang S. and Yao W. Semiparametric mixtures of nonparametric regressions. Annals of the Institute of Statistical Mathematics, 70:131–154, 2018.
2. Benaglia, T., Chauveau, D., Hunter, D. R., & Young, D. S. mixtools: an R package for analyzing mixture models. Journal of statistical software, 32:1-29, 2010.

