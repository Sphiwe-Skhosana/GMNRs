
# Introduction
This folder contains all the ```R``` code used to fit semi-parametric Gaussian mixtures of partially linear models (SPGMPLMs) uisng the proposed objective-based approach, the profile-likelihoo EM (PL-EM) approach proposed by Wu and Liu (2017) and the naive EM algorithm.

# Description of the code
In this section, we provide a brief description of the code (contained in the ```R``` script Fitting-SPGMPLMs.R):
* The ```SPGMPLMs_OB_EM(x,t,y,k,bw,xgrid,init.model)``` function fits the NPGMNRs model using the objective-based approach.
* The ```SPGMPLMs_Naive_EM (x,t,y,k,bw,xgrid,init.model)``` function fits the SPGMPLMs using the naive EM algorithm without consideration to the label-switching problem.
* The ```SPGMPLs_PL_EM(x,t,y,k,bw,xgrid,init.model)``` function fits the SPGMPLMs using the PL_EM algorithm of Wu and Liu (2017).

  #### Arguments (inputs)
  + ```x``` is an $n\times p$ matrix of covariates, excluding the intercept.
  + ```y``` is a vector of length $n$ that consists of the response variable values
  + ```k``` is the number of components
  + ```bw``` is the bandwidth 
  + ```xgrid``` is a vector of length $N\leq n$ that consists of a set of local grid points
  + ```init.model``` is a list object that contains the model to initialize the algorithm. It is obtained as the output of the ```initialize.model``` function. See the description of the ```initialize.model``` function below.
    
  #### Values (outputs)
  All of the above functions return a list with the following items:
  + ```resp``` is an $n\times k$ matrix of the fitted responsibilities
  + ```mix.prop``` is the fitted component mixing proportion functions
  + ```mix.mu``` is the fitted component regression functions (CRFs)
  + ```mix.sigma2``` is the fitted component variance functions
  + ```LL``` is the fitted log-likelihood
* The ```initialize.model(x,y,k,method=1,true.init.model=NULL,p=1)``` function fits/computes a mixture of regressions model using one of two methods specified by the argument ```method```. The fitted model is used to initialize any of the above functions used to fit the NPGMNRs model.
    #### Arguments (inputs)
  + ```x``` is a vector of length $n$ that consists of the covariate values. Note that the function can only take only one covariate.
  + ```y``` is a vector of length $n$ that consists of the response variable values
  + ```k``` is an integer that specifies the number of components
  + ```method``` is an integer that specifies the initialization method that should be used. The choices are ```1 - mixture of polynomial regressions``` and ```2 - the true NPGMNRs model (if known)```
  + ```true.init.model``` is a list object that contains the true model
  + ```p``` If ```method=1``` then this parameter specifies the degree of the polynomial.    
  #### Values (output)
  The ```initialize.model``` function returns a list object ```init.model0``` with the initial fitted model.
# References
1. X, Wu. and T, Liu. Estimation and testing for semiparametric mixtures of partially linear models. Communications in Statistics-Theory and Methods, 46(17), 8690-8705, 2017

