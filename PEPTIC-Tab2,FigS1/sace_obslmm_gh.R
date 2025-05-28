
##data 
#source("/home/User/sacediagreal_final_gh.R")
df_peptic_surv<-df_peptic_trunc_host_sub_ord[df_peptic_trunc_host_sub_ord$fail_hosp_event==1,]
df_peptic_surv$period<-df_peptic_surv$period-1

#lmm model fit on observed
fit_surv<-lme4::lmer(log(HOS_LOS)~trtA*AGE+trtA*GENDER+trtA*AP2_ADMIT+period+(1|clus)+(1|period:clus),data=df_peptic_surv)

coeffs_surv<-fixef(fit_surv)

covs_surv<-cbind(int=1,df_peptic_surv[,c("trtA","AGE","GENDER","AP2_ADMIT","period")])

#SACE diff in means
diff_mean_surv<-coeffs_surv[2]+mean(as.matrix(covs_surv[,c("AGE","GENDER","AP2_ADMIT")])%*%coeffs_surv[7:9])

trt0_mean_exp<-mean(exp(as.matrix(covs_surv[,c("AGE","GENDER","AP2_ADMIT","period")])%*%coeffs_surv[3:6]))

trt1_mean_exp<-mean(exp(coeffs_surv[2]+as.matrix(covs_surv[,c("AGE","GENDER","AP2_ADMIT","period")])%*%coeffs_surv[3:6]+as.matrix(covs_surv[,c("AGE","GENDER","AP2_ADMIT")])%*%coeffs_surv[7:9]))

#SACE relative risk
rr_mean_surv<-trt1_mean_exp/trt0_mean_exp


##cluster bootstrap,
set.seed(01051990)
boot_its<-1000
cl_boot_real<-function(boot){
  idval<-length(unique(df_peptic_surv$clus))
  clustidv<-sample(1:idval,idval,replace=T)
  dflist<-vector("list",idval)
  for(j in 1:idval){
    dflist[[j]]<-df_peptic_surv[df_peptic_surv$clus==clustidv[j],] #resampling by cluster
  }
  dfb<-do.call(rbind,dflist)
  
  fit_survb<-lme4::lmer(log(HOS_LOS)~trtA*AGE+trtA*GENDER+trtA*AP2_ADMIT+period+(1|clus)+(1|period:clus),data=dfb)
  
  coeffs_survb<-fixef(fit_survb)
  
  covs_survb<-cbind(int=1,dfb[,c("trtA","AGE","GENDER","AP2_ADMIT","period")])
  
  diff_mean_survb<-coeffs_survb[2]+mean(as.matrix(covs_survb[,c("AGE","GENDER","AP2_ADMIT")])%*%coeffs_survb[7:9])
  
  #diff_mean_surv_vec[i]<-diff_mean_survb
  
  trt0_mean_expb<-mean(exp(as.matrix(covs_survb[,c("AGE","GENDER","AP2_ADMIT","period")])%*%coeffs_survb[3:6]))
  
  trt1_mean_expb<-mean(exp(coeffs_survb[2]+as.matrix(covs_survb[,c("AGE","GENDER","AP2_ADMIT","period")])%*%coeffs_survb[3:6]+as.matrix(covs_survb[,c("AGE","GENDER","AP2_ADMIT")])%*%coeffs_survb[7:9]))
  
  rr_mean_survb<-trt1_mean_expb/trt0_mean_expb
  
  #rr_mean_surv_vec[i]<-rr_mean_survb
  
  return(c(diff_mean_survb,rr_mean_survb))
  
}

boot_out<-parallel::mclapply(1:boot_its,FUN=cl_boot_real,mc.cores = 5)
boot_out_mat<-do.call(rbind, boot_out)
diff_mean_surv_vec<-boot_out_mat[,1]
quantile(diff_mean_surv_vec,c(0.025,.975))
rr_mean_surv_vec<-boot_out_mat[,2]
quantile(rr_mean_surv_vec,c(0.025,.975))