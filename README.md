# Estimating growth curve parameters via MCMC in JAGS

* Treating parameters as independent, and estimating parameters and their variance: `iid_param.Rmd`
* Treating parameters as covarying, and estimating parameters and their variance/covariance: `covar_param.Rmd` 
  * refer `e_outliers.Rmd` for removal of outlier genotypes to improve estimates for e
* Treating paramaters as independent, using Cholesky factorisation to incorporate genotype relationship matrix (GRM) and estimating parameters: `GRM_iidparam.Rmd`
* Incorporating a Kronecker product of covariance matrix and GRM in the Cholesky factorisation: `GRM_chol.R`
  * Uses `runjags`, with add-on `cholt` function (using [extendJAGS](https://github.com/meerachotai/extendJAGS))
  * TODO: add Kronecker product function from LAPACK (using [extendJAGS](https://github.com/meerachotai/extendJAGS))
