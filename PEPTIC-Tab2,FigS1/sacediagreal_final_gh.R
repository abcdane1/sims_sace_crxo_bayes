require(haven)
require(tidyverse)
require(lme4)
require(coda)

#set random seed
set.seed(18426)
#source("/home/User/crxogibbslogit_finalreal2_gh.R") analysis file 
#peptic file upload
df_peptic<-read_dta("peptic_final_Michael_Harhay_Nov_4_2024.dta")
df_peptic_trunc_host<-df_peptic #renaming to clean up, restrict to key variables

#changes indices 0,1 where 1 is survive until discharge 
df_peptic_trunc_host$fail_hosp_event<-2-df_peptic_trunc_host$fail_hosp_event
#make gender binary, numeric
df_peptic_trunc_host$GENDER<-ifelse(df_peptic_trunc_host$GENDER=="M",1,0)

#variable names
surv_name_rd<-"fail_hosp_event"
out_name_rd<-"HOS_LOS"
log_out_name_rd<-"log(HOS_LOS)"
trt_name_rd<-"trtA"
covs_name_rd<-as.character(substitute(c(AGE,GENDER,AP2_ADMIT)))[-1]
id_name_rd<-"clus"
period_name<-"period"
period_name_rd_cat<-"as.factor(period)" 

df_peptic_trunc_host_sub<-df_peptic_trunc_host[,c(surv_name_rd,out_name_rd,trt_name_rd,covs_name_rd,id_name_rd,"period")]
df_peptic_trunc_host_sub<-df_peptic_trunc_host_sub[complete.cases(df_peptic_trunc_host_sub),]

#ordered nicely
df_peptic_trunc_host_sub_ord<-df_peptic_trunc_host_sub[order(df_peptic_trunc_host_sub$period, df_peptic_trunc_host_sub$clus),]

#fitting initialization models
fixed_effects<-c(trt_name_rd,covs_name_rd,period_name_rd_cat)
random_effects<-c(paste("(","1","|",id_name_rd,")"),paste("(","1","|",id_name_rd,":",period_name,")"))
form_out<-paste(log_out_name_rd,"~",paste(fixed_effects, collapse = " + "), "+", paste(random_effects, collapse = " + "))
form_out<-as.formula(form_out)
form_surv<-paste(surv_name_rd,"~",paste(fixed_effects[-1], collapse = " + "), "+", paste(random_effects, collapse = " + "))
form_surv<-as.formula(form_surv)

#fitted models
#survival status
fit_surv<-lme4::glmer(form_surv,data=df_peptic_trunc_host_sub_ord,family="binomial")
#outcome model 
fit_out<-lme4::lmer(form_out,data=df_peptic_trunc_host_sub_ord)

#variances of res for initialization
#survival
vcov_fit_surv<-as.data.frame(VarCorr(fit_surv))$vcov
#outcome 
vcov_fit_out<-as.data.frame(VarCorr(fit_out))$vcov

#user functions for initial values
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
 
 stratacomp_fun<-function(v){
   p00<-length(which(v==1))/length(v)
   p10<-length(which(v==2))/length(v)
   p11<-1-p00-p10
   return(c(p00,p10,p11))
 }

#initialization parameters (will be uniformly drawn from .5*,1.5*)
sdscz<-sqrt(lognormal_var(vcov_fit_surv[2],0))
sdscpz<-sqrt(lognormal_var(vcov_fit_surv[1],0))
sdscw<-sdscz
sdscpw<-sdscpz
sdyc11<-sqrt(lognormal_var(vcov_fit_out[2],0))
sdycp11<-sqrt(lognormal_var(vcov_fit_out[1],0))
sdyc10<-sdyc11
sdycp10<-sdycp11
sigma2y11<-lognormal_var(vcov_fit_out[3],0)
sigma2y10<-lognormal_var(vcov_fit_out[3],0)

#thinning parameter
thin<-1

#num chains 
chains<-4

#MCMC total and burn in
total_iterations<-10000
burn_in_iterations<-2500

#number of iterations
included_iterations<-ceiling((total_iterations-burn_in_iterations)/thin)

start<-Sys.time()
param_mat<-sace_truecrxologitv3_noop_full_real(data=df_peptic_trunc_host_sub_ord,trt=trt_name_rd,surv=surv_name_rd,out=out_name_rd,ind=covs_name_rd,clustid="clus",clustp="period",cpsize=NULL,sec=T,het=T,total_iterations=total_iterations,burn_in_iterations=burn_in_iterations,chains=chains,thin=thin,credbound=c(.025,.975),post_summ=F,sace_only=F,logscale=T,rr=T,cs=T,cps=F,cy=T,cpy=T,initv=1)
end<-Sys.time()
# time 
end-start
# params of interest 
key_param_mat<-param_mat[,c("saceest","saceRRest","p00","p10","p11","sigma2_up11","sigma2c11_up","sigma2cp11_up")]

#chains for diagnostics
chain_array<-array(NA,c(included_iterations,chains,ncol(key_param_mat)),dimnames = list(NULL, NULL,colnames(key_param_mat)))

for(i in 1:chains){
chain_array[,i,]<-key_param_mat[(1+included_iterations*(i-1)):(included_iterations*i),]
}

## save to R compatible files
# saveRDS(key_param_mat,"key_param_mat.rds")
# saveRDS(chain_array,"chain_array.rds")

## diagnostics 
# diagnostics<-rstan::monitor(chain_array,warmup=0,thin=thin)
# bayesplot::mcmc_trace(chain_array)

## point estimates
#posterior means
# apply(key_param_mat,2,mean)

## interval estimates
# confidence intervals
# coda:: HPDinterval(coda::as.mcmc(key_param_mat), prob = 1-2*credbound[1])