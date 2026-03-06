#required files

#if same variance across mixture components
#source("/home/User/crxodatasimlogitsepv2_mix_gh.R")
#if different variance across mixture components
#source("/home/User/crxodatasimlogitsepv2_mix_shift_gh.R")
#source("/home/User/crxogibbslogit_finalv4_gh.R")

#number of considered metrics
rr<-T
nmetrics<-6+3*rr

#this is to check failures of convergence of MCMC and return NA  
check_fail <- function(func, ...) {
  result <- tryCatch(
    {
      #call function with args
      func(...)
    },
    error = function(e) {
      #return NA in case of an error
      return(rep(NA,nmetrics))
    }
  )
  return(result)
}

#computing SACE  
sacecrxo_logitwrapper<-function(replic,total_iterations,burn_in_iterations,
                                I,J,csl,csu,mux,zp,wp,
                                y111_1,y110_1,y101_1,y111_2,y110_2,y101_2,
                                sigma2y11,sigma2y10,
                                sdscz,sdscpz,sdscw,sdscpw,
                                sdyc11,sdycp11,sdyc10,sdycp10,mix_par,sec,
                                trtp,logscale,rr,initv){
  
  dfs<-crxodfsimlogit_v2_mix(I=I,J=J,csl=csl,csu=csu,mux=mux,zp=zp,wp=wp,
                             y111_1=y111_1,y110_1=y110_1,y101_1=y101_1,
                             y111_2=y111_2,y110_2=y110_2,y101_2=y101_2,
                             sigma2y11=sigma2y11,sigma2y10=sigma2y10,
                             sdscz=sdscz,sdscpz=sdscpz,sdscw=sdscw,sdscpw=sdscpw,
                             sdyc11=sdyc11,sdycp11=sdycp11,sdyc10=sdyc10,sdycp10=sdycp10,mix_par=mix_par,
                             sec=sec,trtp=trtp,logscale=logscale)
  dfsim<-dfs[[1]]
  
  #keep default names for values, 95% confidence intervals
  
  #fits most complex crxo mixed model all re except cp in ps model
  output1<-check_fail(func=sace_truecrxologitv3_noop_full,data=dfsim,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",
                      sec=sec,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,thin=1,chains=1,post_summ=F,sace_only=T,
                      credbound=c(.025,.975),logscale=logscale,rr=rr,cs=T,cps=F,cy=T,cpy=T,initv=initv)
  
  return(c(output1))
}

#solve for lognormal
lognormal_var<-function(sigma2,mu){
  # Calculate variance of the log-normal distribution
  variance<-(exp(sigma2)-1)*exp(2*mu+sigma2)
  return(variance)
}

#solve for sigma2
inv_lognormal<-function(target_var,mu){
  ratio <- target_var/exp(2*mu)
  x<-(1+sqrt(1+4*ratio))/2
  sigma2<-log(x)
  return(sigma2)
}

#icc to variance
icctovar_fun<-function(bpc,wpc,vari){
  if(bpc>wpc){
    stop("BPC cannot be larger than WPC")
  }else{
    varc<--bpc*vari/(wpc-1)
    varcp<-(bpc-wpc)*vari/(wpc-1)
  }
  return(c(varc,varcp))
}

sigma2y11<-rep(4.670774,12) #var on log-scale 1
sigma2y10<-rep(8.692151,12) #var on log-scale 1.25

#BPC.01,WPC.02 and ICCs=.02
var_scen1_ps<-icctovar_fun(.02,.02,pi^2/3)
var_scen1_out11<-icctovar_fun(.01,.02,inv_lognormal(sigma2y11[1],0))
var_scen1_out10<-icctovar_fun(.01,.02,inv_lognormal(sigma2y10[1],0))

var_log_scale_scen1_out11<-lognormal_var(var_scen1_out11,0)
var_log_scale_scen1_out10<-lognormal_var(var_scen1_out10,0)

#BPC .03, WPC .035 and ICC=.035
var_scen2_ps<-icctovar_fun(.035,.035,pi^2/3)
var_scen2_out11<-icctovar_fun(.03,.035,inv_lognormal(sigma2y11[1],0))
var_scen2_out10<-icctovar_fun(.03,.035,inv_lognormal(sigma2y10[1],0))

var_log_scale_scen2_out11<-lognormal_var(var_scen2_out11,0)
var_log_scale_scen2_out10<-lognormal_var(var_scen2_out10,0)

#BPC .05, WPC .1 and ICC=.1
var_scen3_ps<-icctovar_fun(.1,.1,pi^2/3)
var_scen3_out11<-icctovar_fun(.05,.1,inv_lognormal(sigma2y11[1],0))
var_scen3_out10<-icctovar_fun(.05,.1,inv_lognormal(sigma2y10[1],0))

var_log_scale_scen3_out11<-lognormal_var(var_scen3_out11,0)
var_log_scale_scen3_out10<-lognormal_var(var_scen3_out10,0)

#num clusters 
I<-18
#periods
J<-2
#lb cluster-period size
csl<-50
#ub cluster-period size
csu<-150
mux<-c(.75,.25,-.75) #mean of covariates

y111_1<-c(.1,.15,.15,-.5,.7,.05) #coefficients y11_1(1)
y111_2<-c(.4,.35,.35,-.3,.5,.05) #coefficients y11_2(1)

y110_1<-c(.9,.3,-.15,.1,.05) #coefficients y11_1(0)
y110_2<-c(1.2,.5,.05,-.1,.05) #coefficients y11_2(0)

y101_1<-c(.2,.25,-.3,.15,.075) #coefficients y10_1(1)
y101_2<-c(.5,.45,-.1,-.05,.075) #coefficients y10_2(1)

sdscz<-rep(sqrt(c(var_scen1_ps[1],var_scen2_ps[1],var_scen3_ps[1])),4)
sdscpz<-rep(sqrt(c(var_scen1_ps[2],var_scen2_ps[2],var_scen3_ps[2])),4)
sdscw<-sdscz
sdscpw<-sdscpz
sdyc11<-rep(sqrt(c(var_log_scale_scen1_out11[1],var_log_scale_scen2_out11[1],var_log_scale_scen3_out11[1])),4)
sdycp11<-rep(sqrt(c(var_log_scale_scen1_out11[2],var_log_scale_scen2_out11[2],var_log_scale_scen3_out11[2])),4)
sdyc10<-rep(sqrt(c(var_log_scale_scen1_out10[1],var_log_scale_scen2_out10[1],var_log_scale_scen3_out10[1])),4)
sdycp10<-rep(sqrt(c(var_log_scale_scen1_out10[2],var_log_scale_scen2_out10[2],var_log_scale_scen3_out10[2])),4)
initv_vect<-rep(c(1,2,3),4) #initial variance 

zp<-c(.1,.2,-.4,.1,.05) #z logit coefficients 35,25,40
wp<-c(-.1,-.4,-.3,-.1,.025) #w logit coefficients

rr<-T
logscale<-T

#app tables for log=normal mixture mis-specification

#scenarios with mixing parameter varying 
scenariomat<-cbind(csize=I,csl=csl,csu=csu,bpc=rep(c(.01,.03,.05),4),wpc=rep(c(.02,.035,.1),4),mix_par=rep(c(0.9,0.8,0.7,0.6),each=3))


names_results<-c("TrueSACE","BiasMean","EmpSEMean","RMSEMean","CovgSACEMeanHPD",
                 "TrueP00","TrueP10","TrueP11","BiasP00","BiasP10","BiasP11","EmpSEP00","EmpSEP10","EmpSEP11","RMSEP00","RMSEP10","RMSEP11",
                 "TrueSACERR","BiasMeanRR","EmpSEMeanRR","RMSEMeanRR","CovgSACEMeanRRHPD",
                 "ErrorConv","MeanLogY1_11","MeanLogY0_11","MeanY1_11","MeanY0_11")

resultsmat<-matrix(rep(0,nrow(scenariomat)*length(names_results)),nrow=nrow(scenariomat))

colnames(resultsmat)<-names_results

nsim<-1000

total_iterations <- 10000  #total number of iterations including 
burn_in_iterations <- 2500 #burn-in iterations

for(i in (1:nrow(scenariomat))){
  #set seed here (makes reproducible in chunks)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(01051990)
  
  l<-parallel::mcmapply(sacecrxo_logitwrapper,replic=1:nsim,
                        MoreArgs=list(total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,
                                      I=scenariomat[i,"csize"],J=J,csl=scenariomat[i,"csl"],csu=scenariomat[i,"csu"],
                                      mux=mux,zp=zp,wp=wp,y111_1=y111_1,y110_1=y110_1,y101_1=y101_1,
                                      y111_2=y111_2,y110_2=y110_2,y101_2=y101_2,
                                      sigma2y11=sigma2y11[i],sigma2y10=sigma2y10[i],
                                      sdscz=sdscz[i],sdscpz=sdscpz[i],sdscw=sdscw[i],sdscpw=sdscpw[i],
                                      sdyc11=sdyc11[i],sdycp11=sdycp11[i],sdyc10=sdyc10[i],sdycp10=sdycp10[i],mix_par=scenariomat[i,"mix_par"],
                                      sec=T,trtp=.5,logscale=T,rr=rr,initv=initv_vect[i]),
                        mc.cores=60)
  
  
  truedf<-crxodfsimlogit_v2_mix(I=5000,J=2,csl=scenariomat[i,"csl"],csu=scenariomat[i,"csu"],
                                mux=mux,zp=zp,wp=wp,y111_1=y111_1,y110_1=y110_1,y101_1=y101_1,
                                y111_2=y111_2,y110_2=y110_2,y101_2=y101_2,
                                sigma2y11=sigma2y11[i],sigma2y10=sigma2y10[i],
                                sdscz=sdscz[i],sdscpz=sdscpz[i],sdscw=sdscw[i],sdscpw=sdscpw[i],
                                sdyc11=sdyc11[i],sdycp11=sdycp11[i],sdyc10=sdyc10[i],sdycp10=sdycp10[i],mix_par=scenariomat[i,"mix_par"],
                                sec=T,trtp=.5,logscale=T)[[2]]
  
  gtrue<-truedf$g
  gperctrue<-stratacomp_fun(gtrue) #p00,p10,p11 
  
  truedf11<-truedf[truedf$g==3,]
  truesace<-mean(log(truedf11$y1)-log(truedf11$y0))
  truesacerr<-mean(truedf11$y1)/mean(truedf11$y0)
  

  #mean
  sucsaceestim<-which(!is.na(l[1,])) #successful runs 
  nsimmestim<-length(sucsaceestim) #number of successes
  errorconv<-nsim-nsimmestim #number of failures
  saceestim<-l[1,sucsaceestim] #all viable runs
  biasmean<-mean(saceestim)-truesace #bias mean
  rmsemean<-sqrt(mean((saceestim-truesace)^2)) #rmse mean
  mcsemean<-sqrt(1/((nsimmestim-1))*sum((saceestim-mean(saceestim))^2))
  
  #rr, will set rr=T for sims
  if(rr==T){
    saceestimrr<-l[2,sucsaceestim] #all viable runs
    biasrr<-mean(saceestimrr)-truesacerr #bias median 
    rmserr<-sqrt(mean((saceestimrr-truesacerr)^2))
    mcserr<-sqrt(1/((nsimmestim-1))*sum((saceestimrr-mean(saceestimrr))^2))
    covgperc_rr<-sum(l[(2+3*rr),sucsaceestim] <= truesacerr &  l[(3+3*rr),sucsaceestim] >= truesacerr)/nsimmestim
  } #rmse median
  
  #coverage of posterior HDP interval  
  covgperc<-sum(l[(2+rr),sucsaceestim] <= truesace &  l[(3+rr),sucsaceestim] >= truesace)/nsimmestim #percentile
  
  gpercestim<-l[(4+3*rr):nmetrics,sucsaceestim]
  gpercmean<-apply(gpercestim,1,mean)
  biasperc<-gpercmean-gperctrue
  rmseperc<-sqrt(apply((gpercestim-gperctrue)^2,1,mean))
  mcseperc<-sqrt(1/((nsimmestim-1))*apply((gpercestim-gpercmean)^2,1,sum))
  
  results<-c(truesace,biasmean,mcsemean,rmsemean,covgperc, #mean with coverage
             gperctrue,biasperc,mcseperc,rmseperc, #group percentages 
             truesacerr,biasrr,mcserr,rmserr,covgperc_rr,errorconv,
             mean(log(truedf11$y1)),mean(log(truedf11$y0)),mean(truedf11$y1),mean(truedf11$y0)) #true rr

  
  resultsmat[i,]<-results
  
  if(i %% 1 == 0){
    filename<-paste("crxologlogit_sim_table_mix",i,".csv",sep="")
    finalvals<-cbind(scenariomat,resultsmat)
    #print(finalvals)
    write.csv(finalvals,filename)
  }
}

