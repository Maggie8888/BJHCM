
rm(list=ls())
library(readr)
require(foreach)
require(pcaPP)
require(doMC)
require(ipred)
require(survival)
require(copula)
require(prodlim)
require(doParallel)
####################################################################
#####the function is for log likelihood function of the data########
### function argurments: #######
##b1: the parameter for x1 in the hazard funtion for the terminal event; 
##b2: the parameter for x2 in the hazard funtion for the terminal event;
##b3: the parameter for x1 in the hazard funtion for the recurrent event; 
##b4: the parameter for x2 in the hazard funtion for the recurrent event;
##b01, b02: baseline hazards for the terminal event and the recurrent events, respectively; 
##phi: random effect term denoted as w in paper; {W_i}
##theta: correlation parameter for copula function; {gamma_ij}
## Mu: theta_mu average level correlation
##x1, x2: two covariates in the hazard function for the terminal event and the recurrent event;  
##y1, y2: observed time for the terminal event and the recurrent event; 
##s1, s2: failing indicator for the terminal event and the recurrent event; 
##id: subject id;
##N: the sample size of subjects;   
####################################################################
####################################################################
pow=function(x,y){x^y}


LL1=function(b1, b2, b3, b4, b01, b02, phi, theta, x1, x2, y1, y2, s1,s2,id, N){
  #browser()
  lam1=exp(b01+phi+b1*x1[1]+b2*x2[1]);
  U=exp(-lam1*y1[1]);
  u=lam1*exp(-lam1*y1[1]);
  sum=s1[1]*log(u)+(1-s1[1])*log(U)
  for (i in 1:N){
    lam1=exp(b01+phi+b1*x1[i]+b2*x2[i]);
    lam2=exp(b02+phi+b3*x1[i]+b4*x2[i]);
    U=exp(-lam1*y1[i]);
    u=lam1*exp(-lam1*y1[i]);
    V=exp(-lam2*y2[i]);
    v=lam2*exp(-lam2*y2[i]);
    if(s1[i]==1&&s2[i]==1){
      sum=sum+log(pow(U*V,-theta-1))+log(1+theta)+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v);
    }
    else if(s1[i]==0&&s2[i]==0){
      sum=sum-1/theta*log(pow(U,-theta)+pow(V,-theta)-1)-log(U)
    }
    else if(s1[i]==1&&s2[i]==0){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(U);
    }
    else if(s1[i]==0&&s2[i]==1){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)+(-theta-1)*log(V)+log(v)-log(U)
    }
  }
  return(sum)
}

LL2=function(b1, b2, b3, b4, b01, b02, phi,mu,theta,x1, x2, y1, y2, s1,s2,id, N){
  #browser()
  lam1=exp(b01+phi+b1*x1[1]+b2*x2[1]);
  U=exp(-lam1*y1[1]);
  u=lam1*exp(-lam1*y1[1]);
  sum=s1[1]*log(u)+(1-s1[1])*log(U)
  for (i in 1:N){
    
    lam1=exp(b01+phi+b1*x1[i]+b2*x2[i]);
    lam2=exp(b02+phi+b3*x1[i]+b4*x2[i]);
    U=exp(-lam1*y1[i]);
    u=lam1*exp(-lam1*y1[i]);
    V=exp(-lam2*y2[i]);
    v=lam2*exp(-lam2*y2[i]);
    if(s1[i]==1&&s2[i]==1){
      sum=sum+log(pow(U,-mu*exp(theta[i])-1)*pow(V,-mu*exp(theta[i])-1)*(1+mu*exp(theta[i])))+(-1/(mu*exp(theta[i]))-2)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)+log(v);
      
    }
    else if(s1[i]==0&&s2[i]==0){
      sum=sum-1/(mu*exp(theta[i]))*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1);
    }
    else if(s1[i]==1&&s2[i]==0){
      sum=sum+(-1/(mu*exp(theta[i]))-1)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)-(mu*exp(theta[i])+1)*log(U)+log(u);
    }
    else if(s1[i]==0&&s2[i]==1){
      sum=sum+(-1/(mu*exp(theta[i]))-1)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)-(mu*exp(theta[i])+1)*log(V)+log(v)-log(U);
    }
    
  }
  
  return(sum)
}


######################################################
######   generate data and fit the model   ##########
######################################################
################# time-varying  theta ####################


  out<-function(idx,n.sim){
  
  
  n_patient_s=rep(rep(c(100,200,400),each=9),3) ##possible values for sample size, which is set as 100 ,200 or 400##
  theta_s=rep(c(1,1,1,1.5,1.5,1.5,2,2,2),9)              ##possible values for theta, correlation parameter ##########  
  theta=theta_s[idx]              ##for scenario idx=4, theta is 1.5 ##########            
  n_patient=n_patient_s[idx]      ##for scenario idx=4, sample size is 100###  
  J_s=rep(c(1000,1000,1000),27)    ##for scenario idx=4, J is 100###  
  J=J_s[idx]
  UU_s=rep(c(1,0.6,0.4),27)  ## failure rate 
  UU=UU_s[idx] 
  rho_s=rep(c(0.1,0.3,0.5),each=27)
  rho=rho_s[idx]
  mu=theta
  
  b01=0.5  ###beta_0D=0.5######
  b02=1 ###beta_0R=1########
  b1=1             ###beta_D1=1#######
  b2=1             ###beta_D2=1#######
  b3=2             ###beta_R1=2#######
  b4=2             ###beta_R2=2#######
  n_sample=1000   ###the MCMC posterior sample size is 1000####
  n_burn=200     ###the sample size for burn in period is 200####
  thin=10          ###thinning the posterior MCMC sample to lower down the correlation#####
  slice=thin*((n_burn/thin):(n_sample/thin))
  nsim=nsim
  sigma_epsilon=0.1
  
  
  set.seed(n.sim*100)  ###set the seed for simulation ###
  print(n.sim)
  id=rep(1:n_patient)
  phi=rnorm(n_patient,0,0.5)  
  phi=scale(phi,scale=F)  ###random effect w in the paper##
  x1=rbinom(n_patient,1,0.5) ###a contiunous covariate########
  x2=rnorm(n_patient,0,1) ###a binary covariate############
  lam1=exp(b01+phi[id]+x1*b1+x2*b2) ###hazard for the terminal event###
  lam2=exp(b02+phi[id]+x1*b3+x2*b4) ###hazard for the recurrent event###
  u=runif(n_patient)
  t1=-log(u)/lam1    ###generate terminal event### 
  c1=runif(n_patient,0,UU)
  y1=pmin(t1,c1,0.5) ###observed time####
  s1=(t1<=c1) ###failing indicator####
  theta_predict=matrix()
  main=list()
  
  gamma=matrix(NA,nrow=n_patient, ncol=J+1)
  for(n_i in 1:n_patient){
    t2=0
    j=1
    gamma[n_i,1]=rnorm(1,0,sqrt(sigma_epsilon/(1-rho^2)))
    theta=mu*exp(gamma[n_i,1])
    repeat{
      if(j>=2) gamma[n_i,j]=rho*gamma[n_i,j-1]+rnorm(1,0,sqrt(sigma_epsilon))
      
      theta=mu*exp(gamma[n_i,j])
      
      if(y1[n_i]==t1[n_i]) 
      {
        w=runif(1)
        v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
        #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
        r_new=-log(v)/(lam2[n_i])
      }
      
      else if(s1[n_i]==0){
        w=exp(-lam1[n_i]*y1[n_i])
        uv=rCopula(J,claytonCopula(theta))
        v_sample=which(uv[,1]<=w)
        if (length(v_sample)==0 )  uv=rCopula(J,claytonCopula(theta))
        
        v=sample(uv[uv[,1]<=w,2],1,replace=F)
        r_new=-log(v)/(lam2[n_i])
        # w=runif(1)
        # v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
        # #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
        # r_new=-log(v)/(lam2[n_i])
      }
      if(sum(t2)+r_new<=y1[n_i]&length(t2)<J+1) t2=c(t2,r_new)
      else if (sum(t2)+r_new>y1[n_i]&length(t2)<J+1) {
        t2=c(t2,y1[n_i]-sum(t2))
        t2=t2[-1]
        s2=c(rep(1,length(t2)-1),0)
        break
      }
      else if(length(t2)>=J+1){
        t2=t2[-1]
        s2=c(rep(1,length(t2)-1),1)
        break
      }
      j=j+1
      
    }
    s3=1-s2
    s3[length(s2)]=(s1[n_i])
    main[[n_i]]=cbind(id=n_i,y2=t2,s1=s1[n_i],s2,s3,x1=x1[n_i],x2=x2[n_i],y1=y1[n_i],ni=1:length(t2),maxni=length(t2))
  }
  
  main=data.frame(do.call("rbind",main))
  b1_sample=rep(b1,n_sample)
  b2_sample=rep(b2,n_sample)
  b3_sample=rep(b3,n_sample)
  b4_sample=rep(b4,n_sample)
  b01_sample=rep(b01,n_sample)
  b02_sample=rep(b02,n_sample)
  sigma_sample=rep(0.01,n_sample)
  theta_sample=rep(mu,n_sample) 
  
  #dyn.load("/storage/work/m/mjl6416/jointmodeling/01code/copula_joint_model.so")
  dyn.load("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/copula_joint_model.so")
  
 
  scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1) ###step size of MCMC###
  
  
  phiN_est_list=.C("weibul",b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample), 
                   b01=as.double(b01_sample), b02=as.double(b02_sample), 
                   phi=as.double(phi),sigma=as.double(sigma_sample),
                   theta=as.double(theta_sample), x1=as.double(main$x1), x2=as.double(main$x2), 
                   t1=as.double(main$y1),t2=as.double(main$y2), s1=as.integer(main$s1), 
                   s2=as.integer(main$s2), id=as.integer(main$id-1), no_pat_p=as.integer(length(unique(main$id))), n_burn= as.integer(n_burn),  
                   n_sample=as.integer(n_sample),  n=as.integer(length(main$y1)),scale=as.double(scale))
  
  b01_est1=mean(phiN_est_list$b01[slice]) 
  b1_est1=mean(phiN_est_list$b1[slice]) 
  b2_est1=mean(phiN_est_list$b2[slice]) 
  theta_est1=mean(phiN_est_list$theta[slice])
  b3_est1=mean(phiN_est_list$b3[slice])
  b4_est1=mean(phiN_est_list$b4[slice])
  b02_est1=mean(phiN_est_list$b02[slice])
  sigma_est1=mean(phiN_est_list$sigma[slice])
  
  
  
  b1_sample=rep(b1_est1,n_sample)
  b2_sample=rep(b2_est1,n_sample)
  b3_sample=rep(b3_est1,n_sample)
  b4_sample=rep(b4_est1,n_sample)
  b01_sample=rep(b01_est1,n_sample)
  b02_sample=rep(b02_est1,n_sample)
  sigma_sample=rep(sigma_est1,n_sample)
  rho_sample=rep(rho,n_sample)
  theta_sample=rep(na.omit(c(t(gamma))),n_sample)
  sigma_epsilon_sample=rep(sigma_epsilon,n_sample)
  mu_sample=rep(theta_est1,n_sample)
  
  scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1,0.01)
  dyn.load("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/copula_dep.so") ##load c function##
  
  
  result=.C("weibul",b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
            b01=as.double(b01_sample), b02=as.double(b02_sample),
            phi=as.double(phi),sigma=as.double(sigma_sample),mu=as.double(mu_sample),
            theta=as.double(theta_sample),
            rho=as.double(rho_sample), sigma_epsilon=as.double(sigma_epsilon_sample), x1=as.double(main$x1), x2=as.double(main$x2),
            t1=as.double(main$y1),t2=as.double(main$y2), s1=as.integer(main$s1),
            s2=as.integer(main$s2), id=as.integer(main$id-1), no_pat_p=as.integer(length(unique(main$id))), n_burn= as.integer(n_burn),
            n_sample=as.integer(n_sample),  n=as.integer(length(main$y1)),scale=as.double(scale),ni=as.integer(main$ni),maxni=as.integer(main$maxni))
  
  
  b1_est=mean(result$b1[slice])     
  b2_est=mean(result$b2[slice])
  b3_est=mean(result$b3[slice])
  b4_est=mean(result$b4[slice])
  b01_est=mean(result$b01[slice])
  b02_est=mean(result$b02[slice])
  rho_est=mean(result$rho[slice])
  sigma_epsilon_est=mean(result$sigma_epsilon[slice])
  sigma_est=mean(result$sigma[slice])
  mu_est=mean(result$mu[slice])
  
  ######################################################################
  ######   generate new patient data using time-varying model ##########
  ######################################################################
  n_predict=n_patient ###sample size is n_patient###  
  predict_id=1:n_predict
  b01=0.5          ###beta_0D=0.5######
  b02=1            ###beta_0R=1########
  b1=1             ###beta_D1=1#######
  b2=1             ###beta_D2=1#######
  b3=2             ###beta_R1=2#######
  b4=2             ###beta_R2=2#######
  
  
  
  
 
  id=rep(1:n_predict)
  set.seed(n.sim*100)  ###set the seed for simulation ###
  print(n.sim)
  x1=rbinom(n_predict,1,0.5)
  x2=rnorm(n_predict,0,1)
  phi=rnorm(n_predict,0,0.5)
  phi=scale(phi,scale=F)  ###random effect w in the paper##
  lam1_2=exp(b01+phi+x1*b1+x2*b2)
  lam2_2=exp(b02+phi+x1*b3+x2*b4)
  u=runif(n_predict)
  t1_2=-log(u)/lam1_2
  #t2_2=-log(v)/lam2_2
  c1_2=runif(n_predict,0,UU)
  #c2_2=runif(J,0,1)
  y1_2=pmin(t1_2,c1_2,0.5)
  #y2_2=pmin(t2_2,c2_2)
  quantile(y1_2)[2]
  s1_2=(t1_2<=c1_2)
  #s2_2=(t2_2<=c2_2)
  main_predict2=list()
  theta_predict2=matrix()
  sigma_epsilon=0.1
  gamma=matrix(NA,nrow=n_predict, ncol=J+1)
  theta=matrix(NA,nrow=n_predict, ncol=J+1)
  
  for(n_i in 1:n_predict){
    t2_2=0
    t3_2=0
    j=1
    gamma[n_i,1]=rnorm(1,0,sqrt(sigma_epsilon/(1-rho^2)))
    theta[n_i,1]=mu*exp(gamma[n_i,1])
    repeat{
      if(j>=2) 
        gamma[n_i,j]=rho*gamma[n_i,(j-1)]+rnorm(1,0,sqrt(sigma_epsilon))
      theta[n_i,j]=mu*exp(gamma[n_i,j])
      
      if(y1_2[n_i]==t1_2[n_i]) 
      {
        thetas=theta[n_i,j]
        w=runif(1)
        v=((w^(-thetas/(thetas+1))-1)*(u[n_i]^(-thetas))+1)^(-1/thetas)
        #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
        r_new=-log(v)/(lam2_2[n_i])
      }
      
      else if(s1_2[n_i]==0){
        thetas=theta[n_i,j]
        w=exp(-lam1_2[n_i]*y1_2[n_i])
        uv=rCopula(J,claytonCopula(thetas))
        v_sample=which(uv[,1]<=w)
        if (length(v_sample)==0 )  uv=rCopula(J,claytonCopula(thetas))
        
        v=sample(uv[uv[,1]<=w,2],1,replace=F)
        r_new=-log(v)/(lam2_2[n_i])
        # w=runif(1)
        # v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
        # #v=-1/theta*log(1+w*(1-exp(-theta))/(w*(exp(-theta*u[n_i])-1)-exp(-theta*u[n_i])))
        # r_new=-log(v)/(lam2[n_i])
      }
      if(sum(t2_2)+r_new<=y1_2[n_i]&length(t2_2)<J+1) {
        t2_2=c(t2_2,r_new)
        t3_2=c(t3_2,sum(t2_2))
      }
      else if (sum(t2_2)+r_new>y1_2[n_i]&length(t2_2)<J+1) {
        t2_2=c(t2_2,y1_2[n_i]-sum(t2_2))
        t2_2=t2_2[-1]
        s2_2=c(rep(1,length(t2_2)-1),0)
        t3_2=c(t3_2)
        break
      }
      else if(length(t2_2)>=J+1){
        t2_2=t2_2[-1]
        s2_2=c(rep(1,length(t2_2)-1),1)
        t3_2=c(t3_2)
        break
      }
      j=j+1
      
    }
    s3_2=1-s2_2
    s3_2[length(s2_2)]=(s1_2[n_i])
    main_predict2[[n_i]]=cbind(id=n_i,y2_2=t2_2,s1_2=s1_2[n_i],s2_2,s3_2,x1=x1[n_i],x2=x2[n_i],y1_2=y1_2[n_i],ni=1:length(t3_2),maxni=length(t3_2),t3_2=t3_2)
    theta_predict2=gamma
  }
  
  
  main_predict2=data.frame(do.call("rbind",main_predict2))
  theta_sample=rep(na.omit(c(t(theta_predict2))),n_sample)
  ##########################################################################
  ##########begin prediction and plot the curve of the Brier Scores ########
  ##########################################################################
  ##########################################################################
  ##########begin prediction and plot the curve of the Brier Scores ########
  ##########################################################################
  
  
  t_s=c(0.03,0.06,0.09)
  bs_t<-list()
  bs_t[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_terminal<-list()
  bs_t_terminal[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_terminal[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_terminal[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  for(t_prime in t_s){
    
    t_f=seq(t_prime,0.10,0.01)
    idd<-which(t_s==t_prime)
    for (t_final in t_f){
      set.seed(n.sim*100) 
      s1_t=c()
      predict_terminal=c()
      y1_t=c()
      s1_t2=c()
      predict2=c()
      y1_t2=c()
      for(k in 1:n_predict){
        
        
        main_predict_i2=main_predict2[main_predict2$id==k,]
        if(nrow(main_predict_i2)==1|main_predict_i2$y1_2[1]<t_prime) next;
        s1_t2[k]=(t1_2[k]<c1_2[k])
        y1_t2[k]=min(t1_2[k],c1_2[k])
        s1_t[k]=(t1[k]<c1[k])
        y1_t[k]=min(t1[k],c1[k])
        main_predict_p2=main_predict_i2[main_predict_i2$t3_2<t_prime,]
        main_predict_p2$y1_2=t_prime
        main_predict_p2$y2_2[nrow(main_predict_p2)]=t_prime-main_predict_p2$t3_2[nrow(main_predict_p2)]
        main_predict_p2$s2_2[nrow(main_predict_p2)]=0
        main_predict_s2=main_predict_i2[nrow(main_predict_i2),]
      
        
        set.seed(n.sim*100) 
     ### time-varying model prediction ####
        phiN_sample_post=phiN_sample=c()
        phiN_sample[1]=0
        phiN_sample_post[1]=0
        set.seed(12345)
        j=1
        for(m in slice){
          print(m)
          #####MCMC begin######
          for(un_i in 2:10){
            phiN=0.1*rnorm(1,0,1)+phiN_sample[un_i-1]
            #phiN=rnorm(1,0,result$sigma[m])
            new_ll=LL2(result$b1[m], result$b2[m], result$b3[m], result$b4[m], result$b01[m], result$b02[m], phiN, result$mu[m],
                       result$theta[m], main_predict_p2$x1, main_predict_p2$x2, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s3,main_predict_p2$s2,
                       main_predict_p2$id, nrow(main_predict_p2))+dnorm(phiN,0,result$sigma[m])
            old_ll=LL2(result$b1[m],result$b2[m], result$b3[m], result$b4[m], result$b01[m], result$b02[m], phiN_sample[un_i-1], result$mu[m],
                       result$theta[m], main_predict_p2$x1, main_predict_p2$x2, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s3,main_predict_p2$s2,
                       main_predict_p2$id, nrow(main_predict_p2))+dnorm(phiN_sample[un_i-1],0,result$sigma[m])
            u<-runif(1,0,1)
            log_ratio=new_ll-old_ll
            ###new likelihood######
            if(is.na(log_ratio)) 
            { 
              phiN_sample[un_i]=phiN
            }
            else if(log(u)<log_ratio){
              phiN_sample[un_i]=phiN## accept##
            }
            else   phiN_sample[un_i]= phiN_sample[un_i-1] ##reject$$$
          }
          phiN_sample_post[j]=phiN_sample[10] ###keep the final random effect sample into the random effects for thet(m)###
          j=j+1  
        }
        
        
        phiN_est2=mean(phiN_sample_post)
        predict2[k]=exp(-t_final*exp(b1_est*main_predict_p2$x1[1]+b2_est*main_predict_p2$x2[1]+b01_est+phiN_est2))  
        
        ########### time-invarying model prediction ##########

        phiN_sample_post=phiN_sample=c()
        phiN_sample[1]=0
        phiN_sample_post[1]=0
        set.seed(12345)
        j=1
        for(m in slice){
          print(m)
          #####MCMC begin######
          for(un_i in 2:10){
            phiN=0.1*rnorm(1,0,1)+phiN_sample[un_i-1]
            #phiN=rnorm(1,0,result$sigma[m])
            new_ll=LL1(phiN_est_list$b1[m], phiN_est_list$b2[m], phiN_est_list$b3[m], phiN_est_list$b4[m], phiN_est_list$b01[m], phiN_est_list$b02[m], phiN, 
                       phiN_est_list$theta[m], main_predict_p2$x1, main_predict_p2$x2, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s3,main_predict_p2$s2,
                       main_predict_p2$id, nrow(main_predict_p2))+dnorm(phiN,0,phiN_est_list$sigma[m])
            old_ll=LL1(phiN_est_list$b1[m], phiN_est_list$b2[m], phiN_est_list$b3[m], phiN_est_list$b4[m], phiN_est_list$b01[m], phiN_est_list$b02[m], phiN_sample[un_i-1], 
                       phiN_est_list$theta[m], main_predict_p2$x1, main_predict_p2$x2, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s3,main_predict_p2$s2,
                       main_predict_p2$id, nrow(main_predict_p2))+dnorm(phiN_sample[un_i-1],0,phiN_est_list$sigma[m])
            u<-runif(1,0,1)
            log_ratio=new_ll-old_ll
            ###new likelihood######
            if(is.na(log_ratio)) 
            { 
              phiN_sample[un_i]=phiN
            }
            else if(log(u)<log_ratio){
              phiN_sample[un_i]=phiN## accept##
            }
            else   phiN_sample[un_i]= phiN_sample[un_i-1] ##reject$$$
          }
          phiN_sample_post[j]=phiN_sample[10] ###keep the final random effect sample into the random effects for thet(m)###
          j=j+1  
        }
        
        phiN_est=mean(phiN_sample_post)
        #predict2[k]=exp(-t_final*exp(b1_est*main_predict_p2$x1[1]+b2_est*main_predict_p2$x2[1]+b01_est+phiN_est2))  
        predict_terminal[k]=exp(-t_final*exp(b1_est1*main_predict_p2$x1[1]+b2_est1*main_predict_p2$x2[1]+b01_est1+phiN_est))  
      }
      
      predict2=na.omit(predict2)
      censor2=1-na.omit(s1_t2)
      s1_t2=na.omit(s1_t2)
      y1_t2=na.omit(y1_t2)
      
      predict_terminal=na.omit(predict_terminal)
      censor=1-na.omit(s1_t)
      s1_t=na.omit(s1_t)
      y1_t=na.omit(y1_t)
      #predict_time=seq(0.1,max(y1_t),0.01)
      #predict_time2=seq(0.1,max(y1_t2),0.01)
      s_score_terminal=sbrier(Surv(y1_t[1:length(predict_terminal)],s1_t[1:length(predict_terminal)]),as.vector(predict_terminal),btime=t_final)
      s_score=sbrier(Surv(y1_t2,s1_t2),as.vector(predict2),btime =t_final)
      
      bs_t_terminal[[idd]]=rbind(c(t_final,s_score_terminal),bs_t_terminal[[idd]])
      bs_t[[idd]]=rbind(c(t_final,s_score),bs_t[[idd]])
      
      
    }
    
  }
  
  
  return(list(bs_t,bs_t_terminal))
  }


  
#############################

########################################################################

numCores <- 2
registerDoParallel(numCores)
nsim<-2
idx<-19
start_time <- Sys.time()
result<- foreach(n.sim=1:nsim,.errorhandling = "remove") %dopar% {
  out(idx,n.sim)
}
end_time <- Sys.time()
end_time - start_time
setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03preresults")
saveRDS(result,paste0("400_",idx,".rds"))



