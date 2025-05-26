#data = data set
#trt = treatment variable
#surv = survival variable 
#out = outcome variable
#ind = covariates
#clustid = cluster id 
#clustp = cluster period 
#cpsize = cluster period size
#sec = include period fe
#het = include heterogeneous treatment effect among always-survivors
#total_iterations = total MCMC iterations
#burn_in_iterations = total burn-in iterations
#thin = MCMC thinning parameter, if thin=1, no thinning
#chains = MCMC chains
#post_summ = include complete posterior summary of all parameters
#sace_only = include only SACE summary of sace parameters
	#if sace_only and post_summ are false include all raw posterior draws of parameters
#credbound=credible interval bounds
#logscale = outcome on log scale
#rr= return rr if outcome on log scale
#cs,cps:cs-model cluster random effects in ps model, cps-model has cluster and cluster-period random effects in ps model
#cy,cpy:cy-model cluster random effects in outcome model, cpy-model has cluster and cluster-period random effects in outcome model
#initv = an index for initialization of initials 


sace_truecrxologitv3_noop_full<-function(data,trt="a",surv="s",out="y",ind="x",clustid="id",clustp="j",cpsize="nij",sec=F,
                                         het=T,total_iterations=2100,burn_in_iterations=500,thin=1,chains=1,post_summ=F,sace_only=F,
                                         credbound=c(.025,.975),logscale=F,rr=F,cs=T,cps=T,cy=T,cpy=T,
                                         initv=1){
  
  #remove cluster, then remove cluster period
  if(cs==F){
    cps<-F
  }
  
  if(cy==F){
    cpy<-F
  }
  
  #observed data frame 
  dfsim<-data[,c(trt,surv,out,ind,clustid,clustp,cpsize)]
  
  #renaming data frame
  covl<-length(ind)
  xlabs<-paste(rep("x",covl),1:covl,sep="")
  colnames(dfsim)<-c("a","s","y",xlabs,"id","j","nij")
  
  #iterations
  included_iterations<-ceiling((total_iterations-burn_in_iterations)/thin)
  
  #re-calibrate groups, j=2
  periodind<-dfsim$j-1
  dfsim<-data.frame(dfsim,int=1,jind=periodind)
  
  #individuals 
  n<-nrow(dfsim)
  
  #clusters
  I<-length(unique(dfsim$id)) #may not be by number  
  J<-max(dfsim$j)
  
  #cluster periods, two periods
  IM<-I*J
  
  ##required data needed
  Y<-dfsim$y
  
  if(logscale==T){
    Y<-log(Y)
  }
  
  #id vector
  clustid<-dfsim$id
  
  #cluster period vector
  clustp<-dfsim$j
  
  #cluster period index
  clustpind<-clustid+I*(clustp-1)
  
  #cluster period treatment 1 index
  #sort necessary if not ordered correctly
  clustpa1<-sort(unique(clustpind[dfsim$a==1]))
  
  single_chain_fun<-function(ch){
    
    ##survival strata initial values
    
    #initial values for multinomial logistic parameters
    parmult<-length(xlabs)+1+sec
    thetaz_up<-rep(0,parmult) 
    thetaw_up<-rep(0,parmult) 
    
    #initial values of cluster rand int survival
    etaz_up<-rep(0,I) 
    etaw_up<-rep(0,I) 
    
    #initial values of cluster period rand int survival
    nuz_up<-rep(0,IM) 
    nuw_up<-rep(0,IM) 
    
    #note: *must specify for simulation*
    
    #*initial variance of cluster rand int for survival*
    tau2cz_up<-(sdscz^2)[initv] 
    tau2cw_up<-(sdscw^2)[initv] 
    
    #*initial variance of cluster period rand int for survival*
    tau2cpz_up<-(sdscpz^2)[initv] 
    tau2cpw_up<-(sdscpw^2)[initv] 
    
    ##non-mortal outcome initial values
    
    #*initial variance of non-mortal outcome models*
    sigma2_up11<-ifelse(logscale==F,sigma2y11[initv],inv_lognormal(sigma2y11[initv],0)) #updated from 1 
    sigma2_up10<-ifelse(logscale==F,sigma2y10[initv],inv_lognormal(sigma2y10[initv],0))
    
    #init value of outcome cluster period rand int non-mortal
    xi11_up<-rep(0,I) 
    xi10_up<-rep(0,I) 
    
    #init value of outcome cluster period rand int non-mortal, 
    gamma11_up<-rep(0,IM) 
    gamma10_up<-rep(0,IM) 
    
    #*init value cluster rand int variance non-mortal*
    sigma2c11_up<-ifelse(logscale==F,sdyc11^2,inv_lognormal((sdyc11^2)[initv],0))
    
    #*init value cluster period rand int variance non-mortal*
    sigma2cp11_up<-ifelse(logscale==F,sdycp11^2,inv_lognormal((sdycp11^2)[initv],0))
    
    #*init value cluster rand int variance non-mortal*
    sigma2c10_up<-ifelse(logscale==F,sdyc10^2,inv_lognormal((sdyc10^2)[initv],0))
    #*init value cluster period rand int variance non-mortal*
    sigma2cp10_up<-ifelse(logscale==F,sdycp10^2,inv_lognormal((sdycp10^2)[initv],0))
    
    #initial stratum membership
    g_up<-rep(0,n) 
    
    trts<-dfsim$a
    survs<-dfsim$s
    
    g_up[trts==0 & survs==1]<-3 #survive under treatment 0, 11
    g_up[trts==1 & survs==0]<-1 #die under treatment 1, 00
    g_up[trts==1 & survs==1]<-sample(2:3, sum(trts==1 & survs==1), replace = TRUE) #survive under treatment 1, 10 or 11
    g_up[trts==0 & survs==0]<-sample(1:2, sum(trts==0 & survs==0), replace = TRUE) #die under treatment 0, 00 or 10
    
    #binary indicators 
    g11_up<-ifelse(g_up==3,1,0)
    g10_up<-ifelse(g_up==2,1,0)
    
    ##MCMC runs 
    
    #initiating sequence
    k<-1
    
    for (iter in 1:total_iterations){
      
      ##strata 11 data
      id11<-which(g_up==3)
      D11<-as.matrix(dfsim[id11, c("a","int", xlabs, if (sec==T) "jind")])
      #D11<-as.matrix(dfsim[id11,c("a","int",xlabs,"jind")])
      #D11<-as.matrix(dfsim[id11,c("a","int",xlabs)]) #no effects j
      
      if(het==F){
        D11int<-D11 #if no heterogeneous effect, interaction terms 0
      }
      
      if(het==T){
        #include interactions for conditional treatment effect 
        trt011<-ncol(D11)
        #interaction terms (excludes j)
        inter<-D11[,1]*D11[,-c(1,2,length(xlabs)+3)]
        colnames(inter)<-paste(colnames(D11[,-c(1,2,length(xlabs)+3)]),"int")
        D11int<-cbind(D11,inter)
      }
      
      Y11<-Y[id11]
      #cluster id 11
      clustid11<-clustid[id11]
      #cluster period 11
      clustp11<-clustpind[id11]
      xi11_vect_up<-xi11_up[clustid11] 
      gamma11_vect_up<-gamma11_up[clustp11]
      
      ##strata 10 data, under treatment 1 
      
      id10y<-which(g_up==2 & trts==1)
      D10y<-as.matrix(dfsim[id10y, c("int", xlabs, if (sec==T) "jind")])
      #D10y<-as.matrix(dfsim[id10y,c("int",xlabs,"jind")])
      #D10y<-as.matrix(dfsim[id10y,c("int",xlabs)]) #if no period effect 
      Y10<-Y[id10y]
      
      #cluster id10y
      clustid10y<-clustid[id10y]
      #cluster-period id10y
      clustp10y<-clustpind[id10y]
      xi10_vect_up<-xi10_up[clustid10y]
      gamma10_vect_up<-gamma10_up[clustp10y]
      
      #cluster index matrix for stratum 11
      Py11<-matrix(rep(0,length(id11)*I),ncol=I)
      ry11<-1:length(id11)
      Py11[cbind(ry11,clustid11)]<-1
      
      #cluster index matrix for stratum 10 with non-mortal outcome
      Py10<-matrix(rep(0,length(id10y)*I),ncol=I)
      ry10<-1:length(id10y)
      Py10[cbind(ry10,clustid10y)]<-1
      
      #all non-mortal outcome cluster indices
      Py<-rbind(Py11,Py10)
      
      if(cpy==T){
        #cluster-period index matrix for stratum 11
        Ly11<-matrix(rep(0,length(id11)*IM),ncol=IM)
        Ly11[cbind(ry11,clustp11)]<-1
        
        
        #cluster-period index matrix for stratum 10 with non-mortal outcome
        Ly10<-matrix(rep(0,length(id10y)*IM),ncol=IM)
        Ly10[cbind(ry10,clustp10y)]<-1
        Ly10sub<-Ly10[,clustpa1]
        }
      
      ##update fixed effects for non-mortal outcome among those in g11
      sigpart_11<-solve(t(D11int)%*%D11int+sigma2_up11*.001*diag(ncol(D11int))) #diffuse prior with variance I_11/.001
      mu_11<-sigpart_11%*%(t(D11int)%*%(Y11-xi11_vect_up-gamma11_vect_up)+0) #prior with 0 mean
      sigma2_11<-sigma2_up11*sigpart_11
      
      theta11_up<-MASS::mvrnorm(n=1,mu=mu_11,Sigma=sigma2_11)
      
      ##provides SACE estimate in this iteration
      if(het==F){
        #mean difference
        saceest<-theta11_up[1] #if homogeneous, just parameter 
        if(rr==T){
          #risk ratio
          sacerrest<-exp(theta11_up[1]) #if rr true, just exponentiating parameter due to collapsibility
        }
      }
      
      if(het==T){
        #mean difference (means averaged over covariates) 
        a11<-theta11_up[1]+D11[,c(xlabs)]%*%theta11_up[-c(1:trt011)] #no interaction with j
        saceest<-mean(a11)
        
        if(rr==T){
          #risk ratio (means averaged over covariates) 
          rrestnum<-mean(exp(theta11_up[1]+D11[,c(xlabs, if(sec==T) "jind")]%*%(theta11_up[c(3:trt011)])
                             +D11[,c(xlabs)]%*%(theta11_up[-c(1:trt011)])))
          rrestden<-mean(exp(D11[,c(xlabs, if (sec==T) "jind")]%*%(theta11_up[c(3:trt011)])))
          sacerrest<-rrestnum/rrestden
        }} 
      
      ##update strata 11 variance non mortal outcome
      shapeterm_out11<-.001+(length(Y11))/2
      scaleterm_out11<-.001+t((Y11-D11int%*%theta11_up-xi11_vect_up-gamma11_vect_up))%*%(Y11-D11int%*%theta11_up-xi11_vect_up-gamma11_vect_up)/2
      
      sigma2_up11<-1/rgamma(1,shape=shapeterm_out11,scale=1/scaleterm_out11)
      
      ##update fixed effects for non-mortal outcome among those in g10
      sigpart_10<-solve(t(D10y)%*%D10y+sigma2_up10*.001*diag(ncol(D10y))) #diffuse prior with variance I_10/.001
      mu_10<-sigpart_10%*%(t(D10y)%*%(Y10-xi10_vect_up-gamma10_vect_up)+0) #prior with mean 0
      sigma2_10<-sigma2_up10*sigpart_10
      
      theta10_up<-MASS::mvrnorm(n=1,mu=mu_10,Sigma=sigma2_10)
      
      ##update strata 10 variance non-mortal outcome
      shapeterm_out10<-.001+(length(Y10))/2
      scaleterm_out10<-.001+t((Y10-D10y%*%theta10_up-xi10_vect_up-gamma10_vect_up))%*%(Y10-D10y%*%theta10_up-xi10_vect_up-gamma10_vect_up)/2
      sigma2_up10<-1/rgamma(1,shape=shapeterm_out10,scale=1/scaleterm_out10)
      
      ##update cluster parameters of non-mortal outcomes
      
      if(cy==T){
        #intercept xi of 11
        sigpart_xi11<-solve(t(Py11)%*%Py11/sigma2_up11+1/sigma2c11_up*diag(ncol(Py11)))
        mu_xi11<-(1/sigma2_up11)*sigpart_xi11%*%(t(Py11)%*%(Y11-D11int%*%theta11_up-gamma11_vect_up))
        sigma2_xi11<-sigpart_xi11
      
        
        xi11_up<-MASS::mvrnorm(n=1,mu=mu_xi11,Sigma=sigma2_xi11)
        
        #variance of xi of 11 
        shapeterm_xi11v<-.001+I/2
        scaleterm_xi11v<-.001+t(xi11_up)%*%xi11_up/2
        sigma2c11_up<-1/rgamma(1,shape=shapeterm_xi11v,scale=1/scaleterm_xi11v)
        
        xi11_vect_up<-xi11_up[clustid11]
        
        #xi of 10
        sigpart_xi10<-solve(t(Py10)%*%Py10/sigma2_up10+1/sigma2c10_up*diag(ncol(Py10)))
        mu_xi10<-(1/sigma2_up10)*sigpart_xi10%*%(t(Py10)%*%(Y10-D10y%*%theta10_up-gamma10_vect_up))
        sigma2_xi10<-sigpart_xi10
        
        xi10_up<-MASS::mvrnorm(n=1,mu=mu_xi10,Sigma=sigma2_xi10)
        
        #variance of xi of 10
        shapeterm_xi10v<-.001+I/2
        scaleterm_xi10v<-.001+t(xi10_up)%*%xi10_up/2
        sigma2c10_up<-1/rgamma(1,shape=shapeterm_xi10v,scale=1/scaleterm_xi10v)
        
        xi10_vect_up<-xi10_up[clustid10y]
      }
      
      ##update cluster period parameters of non-mortal outcome models
      if(cpy==T){
        #intercept gamma of 11
        sigpart_gamma11<-solve(t(Ly11)%*%Ly11/sigma2_up11+1/sigma2cp11_up*diag(ncol(Ly11)))
        mu_gamma11<-(1/sigma2_up11)*sigpart_gamma11%*%(t(Ly11)%*%(Y11-D11int%*%theta11_up-xi11_vect_up))
        sigma2_gamma11<-sigpart_gamma11
        gamma11_up<-MASS::mvrnorm(n=1,mu=mu_gamma11,Sigma=sigma2_gamma11)
        
        #variance gamma of 11
        shapeterm_gamma11v<-.001+IM/2
        scaleterm_gamma11v<-.001+t(gamma11_up)%*%gamma11_up/2
        sigma2cp11_up<-1/rgamma(1,shape=shapeterm_gamma11v,scale=1/scaleterm_gamma11v)
        
        #gamma of 10
        sigpart_gamma10<-solve(t(Ly10sub)%*%Ly10sub/sigma2_up10+1/sigma2cp10_up*diag(ncol(Ly10sub)))
        mu_gamma10<-(1/sigma2_up10)*sigpart_gamma10%*%(t(Ly10sub)%*%(Y10-D10y%*%theta10_up-xi10_vect_up))
        sigma2_gamma10<-sigpart_gamma10
        gamma10_up_sub<-MASS::mvrnorm(n=1,mu=mu_gamma10,Sigma=sigma2_gamma10)
        gamma10_up[clustpa1]<-gamma10_up_sub

        #variance of gamma of 10
        shapeterm_gamma10v<-.001+IM/4
        scaleterm_gamma10v<-.001+t(gamma10_up_sub)%*%(gamma10_up_sub)/2
        sigma2cp10_up<-1/rgamma(1,shape=shapeterm_gamma10v,scale=1/scaleterm_gamma10v)
      }
      
      #survival strata model variables
      #dmatc<-as.matrix(dfsim[,c("int",xlabs,"jind")])
      dmatc<-as.matrix(dfsim[,c("int", xlabs, if (sec==T) "jind")])
      
      etaz_vect_up<-etaz_up[clustid]
      nuz_vect_up<-nuz_up[clustpind]
      etaw_vect_up<-etaw_up[clustid]
      nuw_vect_up<-nuw_up[clustpind]
      
      ##update latent w11_up, polya-gamma
      tilt11_pg<-dmatc%*%thetaz_up+etaz_vect_up+nuz_vect_up-
        log(1+exp(dmatc%*%thetaw_up+etaw_vect_up+nuw_vect_up))
      
      w11_up<-BayesLogit::rpg(num=n,h=1,z=tilt11_pg)
      
      ##update latent w11_up, polya-gamma
      tilt10_pg<-dmatc%*%thetaw_up+etaw_vect_up+nuw_vect_up-
        log(1+exp(dmatc%*%thetaz_up+etaz_vect_up+nuz_vect_up))
      
      w10_up<-BayesLogit::rpg(num=n,h=1,z=tilt10_pg)
      
      ##for multinomial logistic regression parameters
      
      #diffuse prior 
      sigma2_indef<-diag(.001,ncol(dmatc))
      
      #theta_z for stratum 11 params (assumes mean 0 prior)
      cov_thetaz<-solve(t(dmatc*w11_up)%*%dmatc+sigma2_indef)
      mean_thetaz<-cov_thetaz%*%(t(dmatc)%*%(g11_up-1/2)-
                                   t(dmatc*w11_up)%*%(etaz_vect_up+nuz_vect_up-log(1+exp(dmatc%*%thetaw_up+etaw_vect_up+nuw_vect_up))))  
      thetaz_up<-MASS::mvrnorm(n=1,mu=mean_thetaz,Sigma=cov_thetaz)
      
      #theta_w for stratum 10 params (assumes mean 0 prior)
      cov_thetaw<-solve(t(dmatc*w10_up)%*%dmatc+sigma2_indef)
      mean_thetaw<-cov_thetaw%*%(t(dmatc)%*%(g10_up-1/2)-
                                   t(dmatc*w10_up)%*%(etaw_vect_up+nuw_vect_up-log(1+exp(dmatc%*%thetaz_up+etaz_vect_up+nuz_vect_up))))  
      
      thetaw_up<-MASS::mvrnorm(n=1,mu=mean_thetaw,Sigma=cov_thetaw)
      
      #cluster index matrix all individuals
      P<-matrix(rep(0,n*I),ncol=I)
      r<-1:n
      #cp<-clustid+I*(clustp-1)
      P[cbind(r,clustid)]<-1
      
      #cluster period index for all individuals
      L<-matrix(rep(0,n*IM),ncol=IM)
      L[cbind(r,clustpind)]<-1
      
      ##cluster intercept parameters of strata model
      
      if(cs==T){
        #etaz, cluster random intercept for stratum 11 in multinomial logistic
        cov_etaz<-solve(t(P*w11_up)%*%P+diag(1/tau2cz_up,I))
        mean_etaz<-cov_etaz%*%(-t(P*w11_up)%*%(dmatc%*%thetaz_up+nuz_vect_up
                                               -log(1+exp(dmatc%*%thetaw_up+etaw_vect_up+nuw_vect_up)))+t(P)%*%(g11_up-1/2))
        etaz_up<-MASS::mvrnorm(n=1,mu=mean_etaz,Sigma=cov_etaz)
        etaz_vect_up<-etaz_up[clustid]
        
        #updating variance of cluster random intercept for 11
        shapterm_etazv<-.001+I/2 #assume every cluster has 11 for now
        scaleterm_etazv<-.001+t(etaz_up)%*%etaz_up/2
        tau2cz_up<-1/rgamma(1,shape=shapterm_etazv,scale=1/scaleterm_etazv)
        
        #etaw, cluster random intercept for stratum 10 in multinomial logistic
        cov_etaw<-solve(t(P*w10_up)%*%P+diag(1/tau2cw_up,I))
        mean_etaw<-cov_etaw%*%(-t(P*w10_up)%*%(dmatc%*%thetaw_up+nuw_vect_up
                                               -log(1+exp(dmatc%*%thetaz_up+etaz_vect_up+nuz_vect_up)))+t(P)%*%(g10_up-1/2))
        etaw_up<-MASS::mvrnorm(n=1,mu=mean_etaw,Sigma=cov_etaw)
        etaw_vect_up<-etaw_up[clustid]
        
        #updating variance of cluster random intercept for 10
        shapterm_etawv<-.001+I/2 #assume every cluster has 10 for now
        scaleterm_etawv<-.001+t(etaw_up)%*%etaw_up/2
        tau2cw_up<-1/rgamma(1,shape=shapterm_etawv,scale=1/scaleterm_etawv)}
      
      ##cluster period intercept parameters of strata model
      
      if(cps==T){
        #updating nuz, cluster-period random intercept for stratum 11 in multinomial logistic
        cov_nuz<-solve(t(L*w11_up)%*%L+diag(1/tau2cpz_up,IM))
        mean_nuz<-cov_nuz%*%(-t(L*w11_up)%*%(dmatc%*%thetaz_up+etaz_vect_up-log(1+exp(dmatc%*%thetaw_up
                                                                                      +etaw_vect_up+nuw_vect_up)))+t(L)%*%(g11_up-1/2))
        nuz_up<-MASS::mvrnorm(n=1,mu=mean_nuz,Sigma=cov_nuz)
        nuz_vect_up<-nuz_up[clustpind]
        
        #updating variance of cluster-period random intercept for 11
        shapeterm_nuzv<-.001+IM/2  #assumes every cluster period has 11 
        scaleterm_nuzv<-.001+t(nuz_up)%*%nuz_up/2
        tau2cpz_up<-1/rgamma(1,shape=shapeterm_nuzv,scale=1/scaleterm_nuzv)
        
        #updating nuw, cluster-period random intercept for stratum 10 in multinomial logistic
        cov_nuw<-solve(t(L*w10_up)%*%L+diag(1/tau2cpw_up,IM))
        mean_nuw<-cov_nuw%*%(-t(L*w10_up)%*%(dmatc%*%thetaw_up+etaw_vect_up-log(1+exp(dmatc%*%thetaz_up
                                                                                      +etaz_vect_up+nuz_vect_up)))+t(L)%*%(g10_up-1/2))
        nuw_up<-MASS::mvrnorm(n=1,mu=mean_nuw,Sigma=cov_nuw)
        nuw_vect_up<-nuw_up[clustpind]
        
        #updating variance of cluster-period random intercept for 10
        shapeterm_nuwv<-.001+IM/2  #assumes every cluster period has 10
        scaleterm_nuwv<-.001+t(nuw_up)%*%nuw_up/2
        tau2cpw_up<-1/rgamma(1,shape=shapeterm_nuwv,scale=1/scaleterm_nuwv)}
      
      ##update strata membership for those that survive under treatment 1
      indt1s1<-which(trts==1 & survs==1)
      Yt1s1<-Y[indt1s1]
      
      #cluster t1s1
      clustidt1s1<-clustid[indt1s1]
      
      #cluster-period t1s1
      clustpt1s1<-clustpind[indt1s1]
      
      #xi subsetted to t1s1
      xit1s1_up11<-xi11_up[clustidt1s1]
      xit1s1_up10<-xi10_up[clustidt1s1]
      
      #gamma 11 and 10 subset to t1s1
      gammat1s1_up11<-gamma11_up[clustpt1s1]
      gammat1s1_up10<-gamma10_up[clustpt1s1]
      
      #eta, nu 11 and 10 subset to t1s1
      etat1s1_z_up<-etaz_vect_up[indt1s1]
      nut1s1_z_up<-nuz_vect_up[indt1s1]
      etat1s1_w_up<-etaw_vect_up[indt1s1]
      nut1s1_w_up<-nuw_vect_up[indt1s1]
      
      #parameters subset to t1s1 
      D11t1s1<-as.matrix(dfsim[indt1s1, c("a","int", xlabs, if (sec==T) "jind")])
      #D11t1s1<-as.matrix(dfsim[indt1s1,c("a","int",xlabs,"jind")])
      
      if(het==F){
        D11intt1s1<-D11t1s1 #homogeneous treatment effect 
      }
      if(het==T){
        #heterogeneous treatment effect 
        intert1s1<-D11t1s1[,1]*D11t1s1[,-c(1,2,length(xlabs)+3)] #note, does not include j
        colnames(intert1s1)<-paste(colnames(D11t1s1[,-c(1,2,length(xlabs)+3)]),"int")
        D11intt1s1<-cbind(D11t1s1,intert1s1)
      }
      
      #strata covariates t1s1
      D10t1s1<-D11t1s1[,-1]
      #outcome model covariates t1s1
      D10t1s1y<-D11t1s1[,-1] #remove j from outcome model y10
      #D10t1s1y<-D11t1s1[,-c(1,length(xlabs)+3)]
      
      #outcome models for 11 and 10
      mu11<-D11intt1s1%*%theta11_up+xit1s1_up11+gammat1s1_up11
      f11<-dnorm(Yt1s1,mean=mu11,sd=sqrt(sigma2_up11))
      
      mu10<-D10t1s1y%*%theta10_up+xit1s1_up10+gammat1s1_up10
      f10<-dnorm(Yt1s1,mean=mu10,sd=sqrt(sigma2_up10))
      
      #logit values for 11 and 10 
      mt1s1_11<-exp(D10t1s1%*%thetaz_up+etat1s1_z_up+nut1s1_z_up)
      mt1s1_10<-exp(D10t1s1%*%thetaw_up+etat1s1_w_up+nut1s1_w_up)
      
      #probability always survivor 
      p1<-(f11*mt1s1_11)/(f11*mt1s1_11+f10*mt1s1_10)
      
      #strata membership
      gt1s1_up<-rbinom(length(p1),1,p1)+2
      g_up[indt1s1]<-gt1s1_up
      
      #indicators for membership
      g11_up[indt1s1]<-ifelse(gt1s1_up==3,1,0)
      g10_up[indt1s1]<-ifelse(gt1s1_up==2,1,0)
      
      ##update strata membership for those that die under treatment 0
      indt0s0<-which(trts==0 & survs==0)
      Yt0s0<-Y[indt0s0]
      
      #eta and nu 10 subset to t0s0
      etat0s0_w_up<-etaw_vect_up[indt0s0]
      nut0s0_w_up<-nuw_vect_up[indt0s0]
      
      #strata covariates t0s0
      D10t0s0<-as.matrix(dfsim[indt0s0, c("int", xlabs, if (sec==T) "jind")])
      #D10t0s0<-as.matrix(dfsim[indt0s0,c("int",xlabs,"jind")])
      #D10t0s0<-as.matrix(dfsim[indt0s0,c("int",xlabs,"jind")])
      
      #probability model for 10, and 00 (is 1/exp())
      mt0s0_10<-exp(D10t0s0%*%thetaw_up+etat0s0_w_up+nut0s0_w_up)
      p0<-mt0s0_10/(1+mt0s0_10)
      
      #strata membership
      gt0s0_up<-rbinom(length(p0),1,p0)+1
      g_up[indt0s0]<-gt0s0_up
      
      #indicators for membership
      g10_up[indt0s0]<-ifelse(gt0s0_up==2,1,0)
      
      if(iter==1){
        par_names<-c("saceest","p00","p10","p11",paste(names(theta11_up),"th11"),
                     paste(names(theta10_up),"th10"),
                     paste(names(thetaz_up),"thz"),
                     paste(names(thetaw_up),"thw"),"sigma2_up11","sigma2_up10","sigma2c11_up","sigma2cp11_up",
                     "sigma2c10_up","sigma2cp10_up","tau2cz_up","tau2cpz_up","tau2cw_up","tau2cpw_up")
        
        if(rr==T){
          par_names<-c(par_names[1],"saceRRest",par_names[-1])
          
        }
        parameter_upmat<-matrix(rep(0,length(par_names)*included_iterations),ncol=length(par_names))
        colnames(parameter_upmat)<-par_names
      }
      
      if (iter > burn_in_iterations && (iter-burn_in_iterations-1)%%thin == 0) {
        #minimal info for now
        gperc<-stratacomp_fun(g_up)
        out<-c(saceest,gperc,theta11_up,theta10_up,thetaz_up,thetaw_up,sigma2_up11,sigma2_up10,sigma2c11_up,sigma2cp11_up,
               sigma2c10_up,sigma2cp10_up,tau2cz_up,tau2cpz_up,tau2cw_up,tau2cpw_up)
        
        if(rr==T){
          out<-c(out[1],sacerrest,out[-1])
          
        }
        parameter_upmat[k,]<-out
        k<-k+1
      }}
    
    
    return(parameter_upmat)}
  
  if(chains==1){
    parameter_upmat<-single_chain_fun(1)
  }else{
    parameter_uplist<-parallel::mclapply(1:chains,FUN=single_chain_fun,mc.cores=min(chains,6))
    parameter_upmat<-do.call(rbind,parameter_uplist)
  }
  
  if(post_summ==T){
    sace_only<-F
  }
  
  if(post_summ==T){
    #posterior mean
    posterior_par_output<-apply(parameter_upmat,2,mean)
    #posterior median
    #posterior_par_output<-apply(parameter_upmat,2,median)
    #HPD interval
    cred_int_output<-coda:: HPDinterval(coda::as.mcmc(parameter_upmat), prob = 1-2*credbound[1])
    #percentile interval
    #cred_int_output<-t(apply(parameter_upmat, 2, function(x) quantile(x, probs = c(credbound[1], credbound[2]))))
    output<-list(posterior_par_output,cred_int_output)}
  
  if(sace_only==T & post_summ==F){
    posterior_par_output<-apply(parameter_upmat[,1:(1+rr+3)],2,mean)
    if(rr==T){
      cred_int_output<-coda:: HPDinterval(coda::as.mcmc(parameter_upmat[,1]), prob = 1-2*credbound[1])
      cred_int_outputRR<-coda:: HPDinterval(coda::as.mcmc(parameter_upmat[,2]), prob = 1-2*credbound[1])
      names(cred_int_output)<-c("saceestLB","saceestUB")
      names(cred_int_outputRR)<-c("saceestRRLB","saceestRRUB")
      output<-c(posterior_par_output[1:(1+rr)],cred_int_output,cred_int_outputRR,posterior_par_output[-c(1:(1+rr))])
    }else{
      cred_int_output<-coda:: HPDinterval(coda::as.mcmc(parameter_upmat[,1]), prob = 1-2*credbound[1])
      names(cred_int_output)<-c("saceestLB","saceestUB")
      output<-c(posterior_par_output[1:(1+rr)],cred_int_output,posterior_par_output[-c(1:(1+rr))])}
  }
  
  if(post_summ==F & sace_only==F){
    output<-parameter_upmat
  }
  return(output)
}