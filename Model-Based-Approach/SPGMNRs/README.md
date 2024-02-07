# Introduction
This folder contains all the ```R``` code used to fit semi-parametric mixtures of non-parametric regressions (SPGMNRs) uisng the proposed model-based approach, the local EM approach proposed by Xiang and Yao (2018) and the naive EM algorithm.

# Description of the code
In this section, we provide a brief description of the code (contained in the ```R``` script Fitting-SPGMNRs.R):
* The ```SPGMNR_MB_ECM(x,y,k,bw,xgrid,init.model,lmd)``` function fits the SPGMNRs model using the model-based ECM approach.
* The ```SPGMNR_MB_EM (x,y,k,bw,xgrid,init.model,lmd)``` function fits the SPGMNRs model using the model-based EM approach.
* The ```SPGMNR_Naive_EM (x,y,k,bw,xgrid,init.model)``` function fits the SPGMNRs model using the naive EM algorithm without consideration to the label-switching problem.
* The ```SPGMNR_LEM (x,y,k,bw,xgrid,init.model)``` function fits the SPGMNRs model using the local EM algorithm of Xiang and Yao (2018). The LEM algorithm makes use of the effective EM algorithm to account for the label-switching problem.

  #### Arguments (inputs)
  + ```x``` is the covariate value. The function can only take only one covariate.
  + ```y``` is the response variable values
  + ```k``` is the number of components
  + ```bw``` is the bandwidth 
  + ```xgrid``` is the set of local grid points
  + ```init.model``` is a list object that contains the model to initialize the algorithm. See the description of the ```initialize.model``` function.
  + ```lmd``` is the threshold parameter value
    
  #### Values (outputs)
  All of the above functions return a list with the following items:
  + ```resp``` is an $n\times k$ matrix of the fitted responsibilities
  + ```mix.prop``` is the fitted mixing proportions
  + ```mix.mu``` is the fitted component regression functions (CRFs)
  + ```mix.sigma2``` is the fitted variances
  + ```df``` is the total degrees of freedom for the fitted CRFs
  + ```BIC``` is the value of the fitted Bayesian information criterion
  + ```LL``` is the fitted log-likelihood
* The ```initialize.model(x,y,k,method=1,true.init.model=NULL,p=1)``` function fits/computes a mixture of regressions model using one of three methods specified by the argument ```method```. The fitted model is used to initialize the any of the above functions used to fit the SPGMNRs.
# References
1. Xiang S. and Yao W. Semiparametric mixtures of nonparametric regressions. Annals of the Institute of Statistical Mathematics, 70:131–154, 2018.
