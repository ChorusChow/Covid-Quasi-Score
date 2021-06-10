
QSOEID: hybrid statistical-epidemiological models for prediction in transmission disease with application to SARS-CoV-2
==============================================
  
  
## Outline
1. Description
2. QSOEID (Quasi-Score online estimation for infectious disease) workflow
3. Package requirements
4. Run QSOEID
5. Results and running time
6. How to run QSOEID on your data?
  
## Description
This README is prepared for journal peer review of the "hybrid statistical-epidemiological models for prediction in transmission disease with application to SARS-CoV-2" paper. 

The proposed Quasi-Score online estimator is proposed for estimating the instantaneous reproduction number of an transmission disease that meets the basic assumptions of time-since-infection model with daily incident cases and covariates data. It allow covariates with measurement error to participate in the model while impose no distributional assumptions, and can update estimators online whenever new data are available.

To demonstrate its usage and performance (as in Section 5 of the paper), we use a simulated data (incident cases of a transmission disease) with time length T=120, initiated incident cases I_0=500 and two simulated covariates data to testify its performance on estimating the association of the instantaneous reproduction number of a transmission disease and potiential covariates. The current code version is written for a special case of the model (as described in Section 4 of the paper) with the link function h being the log function, I_t|\mathcal{F}_{t-1} \sim Poisson(R_t \Lambda_t), the model structure admits AR(1) form and allows two covariates inputs (required to be invertible for its sample covariance matrix with days greater than two.)

Further, the code provides an bootstrap inspection on construct a bootstrap confidence interval for the parameters and a bootstrap confidence band for the instantaneous reproduction number with default tunning parameter block length tunningl=45.

QSOEID is a function coded in QSOEID.R file under [Covid-Quasi-Score](https://github.com/ChorusChow/Covid-Quasi-Score). 

## QSOEID workflow (Figure 2)
![](workflow.png)

## Package Requirements
- A database with clear and consistent variable names
- R version: R (>= 4.0.2)
- On Windows: download and install [pracma](https://CRAN.R-project.org/package=pracma), [ggplot2](https://CRAN.R-project.org/package=ggplot2), [EpiEstim](https://CRAN.R-project.org/package=EpiEstim), [nlme](https://CRAN.R-project.org/package=nlme), [dlnm](https://CRAN.R-project.org/package=dlnm), [tsModel](https://CRAN.R-project.org/package=tsModel), [corrplot](https://CRAN.R-project.org/package=corrplot), [mvtnorm](https://CRAN.R-project.org/package=mvtnorm), [stats](https://CRAN.R-project.org/package=stats), [gridExtra](https://CRAN.R-project.org/package=gridExtra), [ggridges](https://CRAN.R-project.org/package=ggridges),
[lubridate](https://CRAN.R-project.org/package=lubridate)

## Run QSOEID example with code

In the example below we aim to testify the performance of proposed quasi-score online estimation method. Two daily covariates was generated independently from a normal distribution with trend term and a uniform distribution with logistic transformation. The first covariate is generated to mimic the temperature of Philadelphia between Mar.1st and May.31st, 2020, while the second is for mimicking the social distancing measured by a percent change in visits to nonessential businesses revealed by daily cell-phone movement within each county. We using a AR(1) structure and assume the fitted model is correctly specified and thus the parameters of interest would be denoted as (\phi_0,\theta_1,\beta_1,\beta_2).

We run the example in local directory.

Step 0: load related R packages and set tunning parameters

```r
## load packages
library(pracma)
library(ggplot2)
library(EpiEstim)
library(nlme)
library(dlnm)
library(tsModel)
library(corrplot)
library(mvtnorm)
library(stats)
library(gridExtra)
library(ggridges)
library(lubridate)

## set working directory
setwd("~/Desktop/QSOEID")

### Set tunning parameter
## tau_0, pre-specified time point where before date tau_0, the MLE will be applied. Default tau_0=5.
tau_0=5
## NoCov, number of covariates we choose
NoCov=2
## T, number of observations/days
T=120
### R[0], the instantaneous reproduction number at time 0.
R_0=3
## I_0, start cases
I_0=500
## rep, number of replications for bootstrap
rep=200
## tunningl, the block length of a fragment of the time series (daily incident cases) used in bootstrap
tunningl=45
## bias_corr_const, the bias correction constant, default is bias_corr_const=1
bias_corr_const=exp(-0.001/2)
``` 

Step 1: Sample data. To avoid redundancy, we generated sample in a seperated R-file with oracle parameter values (\phi_0,\theta_1,\beta_1,\beta_2)=(0.5,0.7,-0.02,-0.125), and the sampled data includes covariates Z, being a (T \times NoCov) matrix, and incident cases I, being a (1 \times T) matrix.

```r
load("sampled_data.rda")
```

Step 1: DLMM Initialization

```r
## setup pda control
control <- list(project_name = 'Length of stay study',
                step = 'initialize',
                sites = c('site1', 'site2', 'site3'),
                heterogeneity = TRUE,
                heterogeneity_effect = 'random',
                model = 'DLM',
                family = 'gaussian',
                outcome = "los",
                variables = c('age', 'sex', 'lab'), 
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )


## specify your working directory, default is the current working dir
mydir <- getwd()   
pda(site_id = 'site1', control = control, dir = mydir)
# you now can see control.json in the working dir


## DO NOT RUN: in actual collaboration, account/password for pda server will be assigned, thus:
# pda(site_id = 'site1', control = control, uri = 'https://pda.one', secret='abc123')
## you can also set your environment variables, and no need to specify them in pda:
# Sys.setenv(PDA_USER = 'site1', PDA_SECRET = 'abc123', PDA_URI = 'https://pda.one')
# pda(site_id = 'site1', control = control)


## site3 communicate its AD: after review, enter "1" to allow tranferring AD 
pda(site_id = 'site3', ipdata = LOS_split[[3]], dir=mydir)
# you now can see site3_initialize.json in the working dir

## site2 communicate its AD: after review, enter "1" to allow tranferring AD   
pda(site_id = 'site2', ipdata = LOS_split[[2]], dir=mydir)
# you now can see site2_initialize.json in the working dir

## site1 communicate its AD: after review, enter "1" to allow tranferring AD   
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)
# you now can see site3_initialize.json in the working dir
# all the AD are ready, control.json is also automatically updated to the next step

``` 

DLMM STEP 2: estimation using AD

```r  
pda(site_id = 'site1', ipdata = LOS_split[[1]], dir=mydir)
# you now can see site1_estimate.json in the working dir

## get the estimated results
config <- getCloudConfig(site_id = 'site1', dir=mydir)
fit.dlmm <- pdaGet(name = 'site1_estimate', config = config)
```


Compare with pooled analysis:
  
  ```r 
## fit LMM using pooled data, assuming random intercepts for each site
fit.pool <- lme4::lmer(los~age+sex+lab+(1|site), REML = F, data=LOS)

# fixed effects (intercept, age, sex, lab) and their sd 
cbind(b.pool = round(summary(fit.pool)$coef[,1], 4),
      b.dlmm = c(fit.dlmm$bhat),      
      sd.pool = round(summary(fit.pool)$coef[,2], 4),  
      sd.dlmm = fit.dlmm$sebhat)  

# variance components (var of random intercepts, and random error)
cbind(vc.pool=round(data.frame(summary(fit.pool)$varcor)$vcov, 4),
      vc.dlmm=round(c(fit.dlmm$Vhat, fit.dlmm$sigmahat^2),4) )

# random intercepts (BLUP) of each sites
cbind(u.pool = round(ranef(fit.pool)$site, 4),
      u.dlmm = c(fit.dlmm$uhat))

```


## Results and Running time

Fixed effect (identical to 2 digits):
  
  |             | b.pool  | b.dlmm  | sd.pool | sd.dlmm|
  |-------------|---------|---------|---------|--------|
  | (Intercept) | 8.1428  | 8.1431  | 0.4926  | 0.4928 |
  | ageold      | 0.9217  | 0.9217  | 0.2840  | 0.2848 |
  | ageyoung    | -0.9825 | -0.9825 | 0.3478  | 0.3487 |
  | sexM        | 0.4533  | 0.4532  | 0.2545  | 0.2552 |
  | lab         | 0.1020  | 0.1020  | 0.0043  | 0.0043 |
  
  Variance components (identical to 2 digits):
  
  |  vc.pool    | vc.dlmm |
  |-------------|---------|
  | 0.4764      | 0.4758  |
  | 16.1371     | 16.2183 |
  
  Random effect (BLUP) (identical to 2 digits):
  
  |        | (Intercept) | u.dlmm  |
  |--------|-------------|---------|
  | site1  | -0.0400     | -0.0403 |
  | site2  | 0.8070      | 0.8063  |
  | site3  | -0.7671     | -0.7661 |
  
  
  The running time in each site in this demo is about 4-5 secs. 


## How to run DLMM on your data?

* Get the data ready, which requires no missing values and clear variable names.
* Define the `control` by specifying the model, data, outcome, variables, site names, and local site.
* Directly run `pda` function as illustrated above.