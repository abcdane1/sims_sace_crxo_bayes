###example execution for SACE in CRXOs 

##required files 

#source("/home/User/crxogibbslogit_finalv4_gh.R")

##inputs for SACE estimation using Bayesian methods for CRXOs 

#data = data set
#trt = treatment variable
#surv = survival variable 
#out = outcome variable
#ind = covariates
#clustid = cluster id 
#clustp = cluster period 
#cpsize = cluster period size
#sec = include period fixed effect
#het = include heterogeneous treatment effect among always-survivors
#total_iterations = total MCMC iterations
#burn_in_iterations = total burn-in iterations
#thin = MCMC thinning parameter, if thin=1, no thinning
#chains = MCMC chains
#post_summ = include complete posterior summary of all parameters
#sace_only = include only SACE summary of SACE parameters
#if sace_only and post_summ are false include all raw posterior draws of parameters
#credbound = credible interval bounds
#logscale = outcome on log scale
#rr = return rr if outcome on log scale
#cs,cps:cs-model cluster random effects in ps model, cps-model has cluster and cluster-period random effects in ps model
#cy,cpy:cy-model cluster random effects in outcome model, cpy-model has cluster and cluster-period random effects in outcome model
#initv = an index for initialization of possible initials for variance/sd parameters, allows for checking of multiple initials (di:change on site)
#if only one specified for each, then set initv=1 as is default 

##running function 

#initializing spread terms, chosen based on data generating mechanism
#NOTE: for one way of deriving initial variance reference terms with real data, see Isenberg et al. 2026 

##non-mortality outcome error variances
#always-survivors 
sigma2y11<-rep(4.670774,3)
#protected-patients
sigma2y10<-rep(8.692151,3)

##sds of random intercepts

#principal strata: always-survivors
sdscz<-0.2591142
sdscpz<-0
#principal strata: protected-patients
sdscw<-sdscz
sdscpw<-sdscpz
#outcome model: always-survivors
sdyc11<-0.1017915
sdycp11<-0.1017915
sdyc10<-0.1140245
sdycp10<-0.1140245

##data (data organized by cluster id sequentially for each period sequentially)
#provided online 
load("/home/User/df_tutorial.RData")

set.seed(5) 

##required functions within SACE estimation function 

#stratum composition 
stratacomp_fun<-function(v){
  p00<-length(which(v==1))/length(v)
  p10<-length(which(v==2))/length(v)
  p11<-1-p00-p10
  return(c(p00,p10,p11))
}

#log scaling for variance 
inv_lognormal<-function(target_var,mu){
  ratio <- target_var/exp(2*mu)
  x<-(1+sqrt(1+4*ratio))/2
  sigma2<-log(x)
  return(sigma2)
}

##running function
sace_truecrxologitv3_noop_full(data=df_tutorial,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",sec=T,
                               het=T,total_iterations=100,burn_in_iterations=50,thin=1,chains=1,post_summ=T,sace_only=F,
                               credbound=c(.025,.975),logscale=T,rr=T,cs=T,cps=F,cy=T,cpy=T,initv=1)



