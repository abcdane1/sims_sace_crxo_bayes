##data simulation for CRXO 

#I = number of clusters
#J = number of periods
#csl,csu = cluster-period bounds, lower csl, upper csu
#mux = covariate means 
#params z: zp (int,x,j)
#params w: wp (int,x,j)
#params y g=11, trt = 1: y111 (trt,int,x,j), 
#params y g=11, trt = 0 : y110 (int,x,j)
	#error var sigma2y11
#params y g=10, trt = 1 : y101 (int,x,j)
	#error var sigma2y10
#sd random effects params
	#z clust,clustp-sdscz,sdscpz, w clust,clusp sdscw,sdscpw
	#y g=11,clust,clustp-sdyc11,sdycp11,y g=10,clust,clustp-sdyc10,sdycp10
#sec - include period
#trtp - treatment sequence probability of assign
#logscale - T=outcome is logged

crxodfsimlogit_v2<-function(I,J,csl,csu,mux,zp,wp,
                         y111,y110,y101,sigma2y11,sigma2y10,
                         sdscz,sdscpz,sdscw,sdscpw,
                         sdyc11,sdycp11,sdyc10,sdycp10,sec=T,
                         trtp=.5,logscale=T){
  
  #cluster size period 1
  cs1<-sample(csl:csu,I,replace=T)
  #cluster size period 2
  cs2<-sample(csl:csu,I,replace=T)
  
  #treatment assignment period 1
  trt1<-rbinom(I,1,p=trtp)
  
  #avoid 0,1 clusters assigned to one treatment at time 1 
  data_test<-function(trt1)
    if (sum(trt1)<2|sum(trt1)>I-1) {
      trt1<-rbinom(I,1,p=trtp)
      data_test(trt1)
    } else {
      return(trt1)
    }
  
  trt1<-data_test(trt1)
  
  #treatment assignment period 2
  trt2<-1-trt1
  
  #cluster size in each period
  n1<-sum(cs1)
  n2<-sum(cs2)
  n<-n1+n2
  
  #empty vectors information 
  
  #covariates x1-x3
  x1<-rep(NA,n)
  x2<-rep(NA,n)
  x3<-rep(NA,n)
  
  #survival outcomes 
  s1<-rep(NA,n) #treatment 1
  s0<-rep(NA,n) #treatment 0
  
  #non-mortal outcomes
  y1<-rep(NA,n) #treatment 1, g11 and g10
  y0<-rep(NA,n) #treatment 0, just g11
  
  id<-c(rep(1:I,cs1),rep(1:I,cs2)) #cluster id
  period<-c(rep(1,n1),rep(2,n2)) #cluster period
  nij<-c(rep(cs1,cs1),rep(cs2,cs2)) #cluster period size 
  
  a<-c(rep(trt1,cs1),rep(trt2,cs2)) #treatment assignment
  s<-rep(NA,n) #observed survival
  y<-rep(NA,n) #observed non-mortal
  
  #augmented data 
  g11<-rep(NA,n) #stratum membership, unknown
  g10<-rep(NA,n) #stratum membership, unknown
  g<-rep(NA,n) #stratum membership, unknown #g=1, 00, g=2, 10, g=3, 11
  
  ##random effect 
  #logcalibration 
  inv_lognormal<-function(target_var, mu){
    # Solve for sigma^2
    ratio <- target_var/exp(2*mu)
    x<-(1+sqrt(1+4*ratio))/2
    sigma2<-log(x)
    return(sigma2)
  }
  
  if(logscale==T){
    #put everything on log scale
    log_sd<-sqrt(inv_lognormal(c(sdscz,sdscpz,sdscw,sdscpw,sdyc11,sdycp11,sdyc10,sdycp10)^2,0))
    sdscz<-log_sd[1];sdscpz<-log_sd[2];sdscw<-log_sd[3];sdscpw<-log_sd[4]
    sdyc11<-log_sd[5];sdycp11<-log_sd[6];sdyc10<-log_sd[7];sdycp10<-log_sd[8]
  }

  #induce within cluster correlation, survival
  Tetai11<-rnorm(I,mean=0,sd=sdscz) #cluster level
  #Teta11<-c(rep(Tetai11,cs1),rep(Tetai11,cs2))
  #cluster period
  Tnuij11<-matrix(c(rnorm(I,0,sd=sdscpz),rnorm(I,0,sd=sdscpz)),ncol=2)
  
  #induce within cluster correlation, survival
  Tetai10<-rnorm(I,mean=0,sd=sdscw) #cluster level
  #Teta10<-c(rep(Tetai10,cs1),rep(Tetai10,cs2))
  #cluster period
  Tnuij10<-matrix(c(rnorm(I,0,sd=sdscpw),rnorm(I,0,sd=sdscpw)),ncol=2) #for now the same 
  
  #induce within cluster correlation, non-mortal
  Txii11<-rnorm(I,mean=0,sd=sdyc11)
  Txii10<-rnorm(I,mean=0,sd=sdyc10)
  #Txi<-c(rep(Txii,cs1),rep(Txii,cs2))
  Tgammaij11<-matrix(c(rnorm(I,0,sd=sdycp11),rnorm(I,0,sd=sdycp11)),ncol=2)
  Tgammaij10<-matrix(c(rnorm(I,0,sd=sdycp10),rnorm(I,0,sd=sdycp10)),ncol=2)
  #Tgamma<-rep(Tgammaij,c(cs1,cs2))
  
  #data generation
  c<-1
  
  csmat<-matrix(c(cs1,cs2),ncol=2) #cluster index by period 
  amat<-matrix(c(trt1,trt2),ncol=2) #treatment by period

  for(j in 1:J){
    for(i in 1:I){
      #covariate means, expressed as vectors
      cij<-csmat[i,j]
      cijrange<-c:(c+cij-1)
      
      mean1<-rep(mux[1],cij)
      mean2<-rep(mux[2],cij)
      mean3<-rep(mux[3],cij) 
      
      #variance matrix of cov vector, independent, all same for now, set to I
      
      #covariates 1-3 generation
      if(logscale==F){
        Sigmacov<-diag(cij) #variance of covariates, hard-coded
        valx1<-MASS::mvrnorm(n=1,mu=mean1,Sigma=3*Sigmacov)
        valx2<-MASS::mvrnorm(n=1,mu=mean2,Sigma=3*Sigmacov)
        valx3<-MASS::mvrnorm(n=1,mu=mean3,Sigma=3*Sigmacov)}
      
      if(logscale==T){
        #hard coded variances
        Sigmacov1<-2.5*inv_lognormal(3,mux[1])*diag(cij) #increase the variance by 1.5
        valx1<-MASS::mvrnorm(n=1,mu=mean1,Sigma=Sigmacov1)
        Sigmacov2<-2*inv_lognormal(3,mux[2])*diag(cij) #increase the variance by 1.5
        valx2<-MASS::mvrnorm(n=1,mu=mean2,Sigma=Sigmacov2)
        Sigmacov3<-inv_lognormal(5,mux[3])*diag(cij)
        valx3<-MASS::mvrnorm(n=1,mu=mean3,Sigma=Sigmacov3)
      }
      
      x1[cijrange]<-valx1
      x2[cijrange]<-valx2
      x3[cijrange]<-valx3
      
      #principal stratum generation 
      Tetaii11<-Tetai11[i]
      Tnuiji11<-Tnuij11[i,j]
      
      Tetaii10<-Tetai10[i]
      Tnuiji10<-Tnuij10[i,j]
      
      if(sec==F){
        zp[5]<-0
        wp[5]<-0
        y111[6]<-0
        y101[5]<-0
        y110[5]<-0
      }
      
      #multinomial model 
      zval11<-exp(zp[1]+zp[2]*valx1+zp[3]*valx2+zp[4]*valx3+zp[5]*(j-1)+Tetaii11+Tnuiji11) #expit
      zval10<-exp(wp[1]+wp[2]*valx1+wp[3]*valx2+wp[4]*valx3+wp[5]*(j-1)+Tetaii10+Tnuiji10)
      
      pval11<-zval11/(1+zval11+zval10)
      pval10<-zval10/(1+zval11+zval10)
      
      pmat<-cbind(pval11,pval10,1-pval11-pval10)
      gvalind<-apply(pmat,1,function(prob) rmultinom(1, 1, prob))
      
      g11val<-gvalind[1,]
      g10val<-gvalind[2,]
      
      g11[cijrange]<-g11val
      g10[cijrange]<-g10val
      
      gval<-ifelse(g11val==1,3,ifelse(g10val==1,2,1))
      g[cijrange]<-gval
      
      #survival status by cluster membership 
      s1val<-ifelse(gval==1,0,1) #survive = 1 under a=1 if 10,11
      s1[cijrange]<-s1val
      
      s0val<-ifelse(gval==3,1,0) #survive = 1 under a=0 if 11
      s0[cijrange]<-s0val
      
      if(logscale==F){
        Sigmay11<-sigma2y11*diag(cij)
        Sigmay10<-sigma2y10*diag(cij)
        } #outcome variance, set to I
      
      if(logscale==T){
        Sigmay11<-inv_lognormal(sigma2y11,0)*diag(cij) #outcome variance under lognormal, variance 3days 
        Sigmay10<-inv_lognormal(sigma2y10,0)*diag(cij) #outcome variance under lognormal, variance 3days 
      }
      
      Txiii11<-Txii11[i]
      Tgammaiji11<-Tgammaij11[i,j]
      
      Txiii10<-Txii10[i]
      Tgammaiji10<-Tgammaij10[i,j]
      
      #outcome y111
      meany1val11<-y111[1]+y111[2]+y111[3]*valx1+y111[4]*valx2+y111[5]*valx3+y111[6]*(j-1)+Txiii11+Tgammaiji11 #mean
      #outcome y101
      meany1val10<-y101[1]+y101[2]*valx1+y101[3]*valx2+y101[4]*valx3+y101[5]*(j-1)+Txiii10+Tgammaiji10 #mean
      #outcome y110
      meany0val11<-y110[1]+y110[2]*valx1+y110[3]*valx2+y110[4]*valx3+y110[5]*(j-1)+Txiii11+Tgammaiji11 #mean
      
      y1val11<-MASS::mvrnorm(n=1,mu=meany1val11,Sigma=Sigmay11) #random draw 111
      y1val10<-MASS::mvrnorm(n=1,mu=meany1val10,Sigma=Sigmay10) #random draw 101
      y0val11<-MASS::mvrnorm(n=1,mu=meany0val11,Sigma=Sigmay11) #random draw 110
      
      if(logscale==T){
        y1val11<-exp(y1val11)
        y1val10<-exp(y1val10)
        y0val11<-exp(y0val11)
      }
      
      #y under treatment 1 according to group, truncated if 00
      y1val<-ifelse(gval==1,NA,ifelse(gval==2,y1val10,y1val11))
      
      #input y1 data
      y1[cijrange]<-y1val
      
      #y under treatment 0 according to group, truncated if 10,00
      y0val<-ifelse(gval==3,y0val11,NA)
      
      #input y0 data
      y0[cijrange]<-y0val
      
      #observed data according to randomization
      aij<-amat[i,j]
      
      if(aij==1){
        s[cijrange]<-s1val
        y[cijrange]<-y1val}
      else{
        s[cijrange]<-s0val
        y[cijrange]<-y0val}
      
      c<-c+cij
    }
  }
  
  dfsim<-data.frame(id=id,j=period,nij=nij,x1=x1,x2=x2,x3=x3,a=a,s=s,y=y)
  dftrue<-data.frame(id=id,j=period,nij=nij,x1=x1,x2=x2,x3=x3,s1=s1,s0=s0,y1=y1,y0=y0,g=g,g11=g11,g10=g10)
  
  l<-list(dfsim,dftrue)
  return(l)}

#stratum composition function
stratacomp_fun<-function(v){
  p00<-length(which(v==1))/length(v)
  p10<-length(which(v==2))/length(v)
  p11<-1-p00-p10
  return(c(p00,p10,p11))
}
