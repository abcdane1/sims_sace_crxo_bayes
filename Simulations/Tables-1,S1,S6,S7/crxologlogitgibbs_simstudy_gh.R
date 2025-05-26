#required files

#source("/home/User/crxodatasimlogitsepv2_gh.R")
#source("/home/User/crxogibbslogit_finalv4_gh.R")

#number of considered metrics
rr<-T
nmetrics<-6+3*rr

#set seed 
RNGkind("L'Ecuyer-CMRG")
set.seed(01051990)

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
                          y111,y110,y101,sigma2y11,sigma2y10,
                          sdscz,sdscpz,sdscw,sdscpw,
                          sdyc11,sdycp11,sdyc10,sdycp10,sec,
                          trtp,logscale,rr,initv){
  
  dfs<-crxodfsimlogit_v2(I=I,J=J,csl=csl,csu=csu,mux=mux,zp=zp,wp=wp,
                                   y111=y111,y110=y110,y101=y101,sigma2y11=sigma2y11,sigma2y10=sigma2y10,
                                   sdscz=sdscz,sdscpz=sdscpz,sdscw=sdscw,sdscpw=sdscpw,
                                   sdyc11=sdyc11,sdycp11=sdycp11,sdyc10=sdyc10,sdycp10=sdycp10,sec=sec,
                                   trtp=trtp,logscale=logscale)
  dfsim<-dfs[[1]]
  
  #keep default names for values, 95% confidence intervals
  
  #fits most complex crxo mixed model all re except cp in ps model
  output1<-check_fail(func=sace_truecrxologitv3_noop_full,data=dfsim,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",
                      sec=sec,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,thin=1,chains=1,post_summ=F,sace_only=T,
                      credbound=c(.025,.975),logscale=logscale,rr=rr,cs=T,cps=F,cy=T,cpy=T,initv=initv)
  
  #fits crxo mixed model without cluster random effects in ps model
  output2<-check_fail(func=sace_truecrxologitv3_noop_full,data=dfsim,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",
                      sec=sec,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,thin=1,chains=1,post_summ=F,sace_only=T,
                      credbound=c(.025,.975),logscale=logscale,rr=rr,cs=F,cps=F,cy=T,cpy=T,initv=initv)
  
  #fits crxo mixed model without cluster random effects in outcome model
  output3<-check_fail(func=sace_truecrxologitv3_noop_full,data=dfsim,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",
                      sec=sec,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,thin=1,chains=1,post_summ=F,sace_only=T,
                      credbound=c(.025,.975),logscale=logscale,rr=rr,cs=T,cps=F,cy=T,cpy=F,initv=initv)

  #fits crxo mixed model without cluster random effects in ps model nor cp re in outcome model
  output4<-check_fail(func=sace_truecrxologitv3_noop_full,data=dfsim,trt="a",surv="s",out="y",ind=c("x1","x2","x3"),clustid="id",clustp="j",cpsize="nij",
                      sec=sec,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,thin=1,chains=1,post_summ=F,sace_only=T,
                      credbound=c(.025,.975),logscale=logscale,rr=rr,cs=F,cps=F,cy=T,cpy=F,initv=initv)

  return(c(output1,output2,output3,output4))
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

##scenarios 

sigma2y11<-rep(4.670774,3) #var on log-scale 1
sigma2y10<-rep(8.692151,3) #var on log-scale 1.25

#BPC .01,WPC.02 and ICCs=.02
var_scen1_ps<-icctovar_fun(.02,.02,pi^2/3)
var_scen1_out11<-icctovar_fun(.01,.02,inv_lognormal(sigma2y11[1],0))
var_scen1_out10<-icctovar_fun(.01,.02,inv_lognormal(sigma2y10[1],0))

var_log_scale_scen1_out11<-lognormal_var(var_scen1_out11,0)
var_log_scale_scen1_out10<-lognormal_var(var_scen1_out10,0)

#BPC .03, WPC .035 and ICCs=.035
var_scen2_ps<-icctovar_fun(.035,.035,pi^2/3)
var_scen2_out11<-icctovar_fun(.03,.035,inv_lognormal(sigma2y11[1],0))
var_scen2_out10<-icctovar_fun(.03,.035,inv_lognormal(sigma2y10[1],0))

var_log_scale_scen2_out11<-lognormal_var(var_scen2_out11,0)
var_log_scale_scen2_out10<-lognormal_var(var_scen2_out10,0)

#BPC .05, WPC .1 and ICCs=.1
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
y111<-c(.1,.15,.15,-.5,.7,.05) #coefficients y11(1)
y110<-c(.9,.3,-.15,.1,.05) #coefficients y11(0)
y101<-c(.2,.25,-.3,.15,.075) #coefficient y10(1)

sdscz<-sqrt(c(var_scen1_ps[1],var_scen2_ps[1],var_scen3_ps[1]))
sdscpz<-sqrt(c(var_scen1_ps[2],var_scen2_ps[2],var_scen3_ps[2]))
sdscw<-sdscz
sdscpw<-sdscpz
sdyc11<-sqrt(c(var_log_scale_scen1_out11[1],var_log_scale_scen2_out11[1],var_log_scale_scen3_out11[1]))
sdycp11<-sqrt(c(var_log_scale_scen1_out11[2],var_log_scale_scen2_out11[2],var_log_scale_scen3_out11[2]))
sdyc10<-sqrt(c(var_log_scale_scen1_out10[1],var_log_scale_scen2_out10[1],var_log_scale_scen3_out10[1]))
sdycp10<-sqrt(c(var_log_scale_scen1_out10[2],var_log_scale_scen2_out10[2],var_log_scale_scen3_out10[2]))
initv_vect<-c(1,2,3) #initial variance 

#table main, app 1
zp<-c(.1,.2,-.4,.1,.05) #z logit coefficients 35,25,40
wp<-c(-.1,-.4,-.3,-.1,.025) #w logit coefficients

#table app 25,25,50 S6,S7
# zp<-c(.7,.2,-.35,.1,.05) #z logit coefficients
# wp<-c(.1,-.225,-.1,-.1,.025) #w logit coefficients

rr<-T
logscale<-T

crxo_opt<-c("CRXO_RE_Full","CRXO_RE_Out","CRXO_RE_PS","CRXO_RE_None")

#scenariomat<-expand.grid(modeltype=crxo_opt,csize=I,csl=csl,csu=csu,bpc=0.05,wpc=0.1)
scenariomat<-cbind(expand.grid(modeltype=crxo_opt,csize=I,csl=csl,csu=csu,bpc=c(.01,.03,.05)),wpc=rep(c(.02,.035,.1),each=4))

#scenariomat<-cbind(rbind(scenariomat,scenariomat),csl=rep(csl,each=nrow(scenariomat)),csu=rep(csu,each=nrow(scenariomat)))

names_results<-c("TrueSACE","BiasMean","EmpSEMean","RMSEMean","CovgSACEMeanHPD",
                 "TrueP00","TrueP10","TrueP11","BiasP00","BiasP10","BiasP11","EmpSEP00","EmpSEP10","EmpSEP11","RMSEP00","RMSEP10","RMSEP11",
                 "TrueSACERR","BiasMeanRR","EmpSEMeanRR","RMSEMeanRR","CovgSACEMeanRRHPD",
                 "ErrorConv")

resultsmat<-matrix(rep(0,nrow(scenariomat)*length(names_results)),nrow=nrow(scenariomat))

colnames(resultsmat)<-names_results

nsim<-1000

total_iterations <- 10000  #total number of iterations including 
burn_in_iterations <- 2500 #burn-in iterations

for(i in c(1,2,3)){
  l<-parallel::mcmapply(sacecrxo_logitwrapper,replic=1:nsim,
                        MoreArgs=list(total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,
                                      I=scenariomat[4*i,"csize"],J=J,csl=scenariomat[4*i,"csl"],csu=scenariomat[4*i,"csu"],
                                      mux=mux,zp=zp,wp=wp,y111=y111,y110=y110,y101=y101,
                                      sigma2y11=sigma2y11[i],sigma2y10=sigma2y10[i],
                                      sdscz=sdscz[i],sdscpz=sdscpz[i],sdscw=sdscw[i],sdscpw=sdscpw[i],
                                      sdyc11=sdyc11[i],sdycp11=sdycp11[i],sdyc10=sdyc10[i],sdycp10=sdycp10[i],
                                      sec=T,trtp=.5,logscale=T,rr=rr,initv=initv_vect[i]),
                        mc.cores=60)
  
  truedf<-crxodfsimlogit_v2(I=5000,J=2,csl=scenariomat[4*i,"csl"],csu=scenariomat[4*i,"csu"],
                            mux=mux,zp=zp,wp=wp,y111=y111,y110=y110,y101=y101,
                            sigma2y11=sigma2y11[i],sigma2y10=sigma2y10[i],
                            sdscz=sdscz[i],sdscpz=sdscpz[i],sdscw=sdscw[i],sdscpw=sdscpw[i],
                            sdyc11=sdyc11[i],sdycp11=sdycp11[i],sdyc10=sdyc10[i],sdycp10=sdycp10[i],
                            sec=T,trtp=.5,logscale=T)[[2]]
    
  gtrue<-truedf$g
  gperctrue<-stratacomp_fun(gtrue)
  
  truedf11<-truedf[truedf$g==3,]
  truesace<-mean(log(truedf11$y1)-log(truedf11$y0))
  truesacerr<-mean(truedf11$y1)/mean(truedf11$y0)
  
  for(j in 1:4){
    #mean
    #truesace<-truesacef[j]
    sucsaceestim<-which(!is.na(l[1+nmetrics*(j-1),])) #successful runs 
    nsimmestim<-length(sucsaceestim) #number of successes
    errorconv<-nsim-nsimmestim #number of failures
    saceestim<-l[1+nmetrics*(j-1),sucsaceestim] #all viable runs
    biasmean<-mean(saceestim)-truesace #bias mean
    rmsemean<-sqrt(mean((saceestim-truesace)^2)) #rmse mean
    mcsemean<-sqrt(1/((nsimmestim-1))*sum((saceestim-mean(saceestim))^2))
    
    #rr, will set rr=T for now
    if(rr==T){
    saceestimrr<-l[2+nmetrics*(j-1),sucsaceestim] #all viable runs
    biasrr<-mean(saceestimrr)-truesacerr #bias median 
    rmserr<-sqrt(mean((saceestimrr-truesacerr)^2))
    mcserr<-sqrt(1/((nsimmestim-1))*sum((saceestimrr-mean(saceestimrr))^2))
    covgperc_rr<-sum(l[(2+3*rr)+nmetrics*(j-1),sucsaceestim] <= truesacerr &  l[(3+3*rr)+nmetrics*(j-1),sucsaceestim] >= truesacerr)/nsimmestim
    } #rmse median
    
    #coverage of posterior HDP interval  
    covgperc<-sum(l[(2+rr)+nmetrics*(j-1),sucsaceestim] <= truesace &  l[(3+rr)+nmetrics*(j-1),sucsaceestim] >= truesace)/nsimmestim #percentile
    
    gpercestim<-l[(4+3*rr):nmetrics+nmetrics*(j-1),sucsaceestim]
    gpercmean<-apply(gpercestim,1,mean)
    biasperc<-gpercmean-gperctrue
    rmseperc<-sqrt(apply((gpercestim-gperctrue)^2,1,mean))
    mcseperc<-sqrt(1/((nsimmestim-1))*apply((gpercestim-gpercmean)^2,1,sum))
    
    results<-c(truesace,biasmean,mcsemean,rmsemean,covgperc, #mean with coverage
               gperctrue,biasperc,mcseperc,rmseperc, #group percentages 
               truesacerr,biasrr,mcserr,rmserr,covgperc_rr,errorconv) #true rr
    
    resultsmat[(4*i-4)+j,]<-results
    }
    
  if(i %% 1 == 0){
    filename<-paste("crxologlogit_sim_table1new",i,".csv",sep="")
    finalvals<-cbind(scenariomat,resultsmat)
    #print(finalvals)
    write.csv(finalvals,filename)
  }
}

