rm(list=ls())

library(readr)
require(foreach)
require(pcaPP)
require(doMC)
require(ipred)
require(survival)
require(copula)
require(prodlim)
library(Rcpp)
library(RcppArmadillo)
library(frailtypack)
require(doParallel)
require(Copula.surv)
setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/")
#source("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/BS.fun.R")
sourceCpp('/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/copula_dep_vec.cpp')
sourceCpp('/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/copula_joint_model_vec.cpp')
Stroke<-read_csv("02data/CHS_data/real_st.csv")
Stroke<-as.data.frame(Stroke)
length(unique(Stroke$id)) # 5465
summary(Stroke$x1)
table(Stroke$x2,exclude = NULL)
table(Stroke$x3,exclude = NULL)
table(Stroke$x4,exclude = NULL)
table(Stroke$x5,exclude = NULL)
Stroke$t3<-Stroke$y2
X<-unique(Stroke$id)
for (i in 1:length(X)){
  temp<-Stroke[Stroke$id==X[i],]
  temp1<-temp$y2
  t3<-vector()
  t3[1]<-temp1[1]
  if (length(temp1)>=2){
    for (j in 2:length(temp1)){
      t3[j]<-sum(temp1[1:j])
    } 
  }
  
  Stroke[Stroke$id==X[i],"t3"]<-t3
}


#### check for the correlation time-varying or not ###
work.data<-Stroke[order(Stroke$id),]
cal_tau<-function(j){
  y2<-vector()
  y1<-vector()
  s1<-vector()
  s2<-vector()
  X=unique(work.data$id)
  for (i in 1:length(X)){
    tmp<-work.data[work.data$id==X[i],]
    if (!is.na(tmp$y2[j]))
    {y2[i]<-tmp$y2[j]
    y1[i]<-tmp$y1[j]
    s2[i]<-tmp$s2[j]
    s1[i]<-tmp$s1[j]
    }
  }
  
  y1<-na.omit(y1)
  y2<-na.omit(y2)
  s1<-na.omit(s1)
  s2<-na.omit(s2)
  result<-U2.Clayton(y1,y2,s1,s2)
  return(result[2])
}

tau<-vector()
for (j in 1:5){
  tau[j]<-cal_tau(j)
}

cat(capture.output(print(tau), file="/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03results/output_CHS_Stroke_tau.txt"))

plot(tau, main="Stroke_CHS") 

out<-function(st_final,main_predict2,n.sim){
  id_idx=lapply(unique(st_final$id), function(x){
    idxs=which(st_final$id==x)
    c(min(idxs),max(idxs))
  })
  id_idx_t=do.call("rbind",id_idx)
  j_start=id_idx_t[,1]
  j_end=id_idx_t[,2]
  
  X<-unique(st_final$id)
  for (i in 1:length(X)){
    temp<-st_final[st_final$id==X[i],]
    temp1<-temp$y1
    st_final[st_final$id==X[i],"maxni"]<-rep(length(temp1),length(temp1))
  }
  
  n_sample=2000
  n_burn=1000
  thin=10
  slice=thin*((n_burn/thin):(n_sample/thin))
  
  J=max(table(st_final$id))
  n_patient=length(unique(st_final$id))
  maxni=J
  
  #x1=st_final$mi_st
  #x2=st_final$GEND01.x
  #x3=st_final$RACE01.x
  #x4=st_final$AVZMSYS
  
  
  
  b1_sample=rep(0.10,n_sample)
  b2_sample=rep(0.43,n_sample)
  b3_sample=rep(0.37,n_sample)
  b4_sample=rep(0,n_sample)
  b5_sample=rep(0.14,n_sample)
  b6_sample=rep(0.46,n_sample)
  b7_sample=rep(0.31,n_sample)
  b8_sample=rep(0,n_sample)
  phi=rep(0,n_patient)
  
  b01_sample=rep(-9,n_sample)
  b02_sample=rep(-13,n_sample)
  sigma_sample=rep(20,n_sample)
  theta_sample=rep(1,n_sample)
  scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1)
  
  set.seed(1234)
  
  dyn.load("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/real_cop.so")
  
  result_cop=.C("weibul",b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample), b5=as.double(b5_sample), b6=as.double(b6_sample),b7=as.double(b7_sample),b8=as.double(b8_sample),
                b01=as.double(b01_sample), b02=as.double(b02_sample), 
                phi=as.double(rep(0,length(unique(st_final$id)))),sigma=as.double(sigma_sample),
                theta=as.double(theta_sample), x1=as.double(scale(st_final$x1)), x2=as.double(st_final$x2),
                x3=as.double(st_final$x3), x4=as.double(st_final$x4), 
                t1=as.double(st_final$y1),t2=as.double(st_final$y2), s1=as.integer(st_final$s3),
                s2=as.integer(st_final$s2), id=as.integer(st_final$id-1), no_pat_p=as.integer(length(unique(st_final$id))), n_burn= as.integer(n_burn),  
                n_sample=as.integer(n_sample),  n=as.integer(nrow(st_final)),scale=as.double(scale),j_start=as.integer(j_start),j_end=as.integer(j_end),ni=as.integer(st_final$ni),maxni=as.integer(st_final$maxni))
  
  
  b01_est1=mean(result_cop$b01) 
  b1_est1=mean(result_cop$b1) 
  b2_est1=mean(result_cop$b2)
  b3_est1=mean(result_cop$b3) 
  b4_est1=mean(result_cop$b4)
  b5_est1=mean(result_cop$b5) 
  b6_est1=mean(result_cop$b6)
  b7_est1=mean(result_cop$b7) 
  b8_est1=mean(result_cop$b8)
  b02_est1=mean(result_cop$b02) 
  sigma_est1=mean(result_cop$sigma) 
  mu_est1=mean(result_cop$theta) 
  
  b1_sample=rep(0.10,n_sample)
  b2_sample=rep(0.43,n_sample)
  b3_sample=rep(0.37,n_sample)
  b4_sample=rep(1,n_sample)
  b5_sample=rep(0.14,n_sample)
  b6_sample=rep(0.46,n_sample)
  b7_sample=rep(0.31,n_sample)
  b8_sample=rep(1,n_sample)
  mu=rep(1,n_sample)
  
  b01_sample=rep(-9,n_sample)
  b02_sample=rep(-13,n_sample)
  sigma_sample=rep(0.01,n_sample)
  gammas=matrix(0,nrow=n_patient, ncol=J)
  maxni=J
  theta_sample=na.omit(c(t(gammas)))
  rho_sample=rep(0,n_sample)
  sigma_epsilon_sample=rep(1,n_sample)
  scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1,0.01)
  
  
  dyn.load("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/01code/real_new.so")
  set.seed(1234)
  result_dep=.C("weibul",b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
                b5=as.double(b5_sample), b6=as.double(b6_sample), b7=as.double(b7_sample), b8=as.double(b8_sample), b01=as.double(b01_sample), b02=as.double(b02_sample),
                phi=as.double(phi),sigma=as.double(sigma_sample),mu=rep(mu,n_sample),
                theta=as.double(theta_sample),
                rho=as.double(rho_sample), sigma_epsilon=as.double(sigma_epsilon_sample), x1=as.double(scale(st_final$x1)), x2=as.double(st_final$x2),
                x3=as.double(st_final$x3), x4=as.double(st_final$x4), t1=as.double(st_final$y1),t2=as.double(st_final$y2), s1=as.integer(st_final$s3),
                s2=as.integer(st_final$s2), id=as.integer(st_final$id-1), no_pat_p=as.integer(length(unique(st_final$id))), n_burn= as.integer(n_burn),
                n_sample=as.integer(n_sample),  n=as.integer(length(st_final$y1)),scale=as.double(scale),ni=as.integer(st_final$ni),maxni=as.integer(st_final$maxni),j_start=as.integer(j_start),j_end=as.integer(j_end))
  
  
  b1_est=mean(result_dep$b1)     
  b2_est=mean(result_dep$b2)
  b3_est=mean(result_dep$b3)
  b4_est=mean(result_dep$b4)
  b5_est=mean(result_dep$b5)     
  b6_est=mean(result_dep$b6)
  b7_est=mean(result_dep$b7)
  b8_est=mean(result_dep$b8)
  b01_est=mean(result_dep$b01)
  b02_est=mean(result_dep$b02)
  rho_est=mean(result_dep$rho)
  sigma_epsilon_est=mean(result_dep$sigma_epsilon)
  sigma_est=mean(result_dep$sigma)
  mu_est=mean(result_dep$mu)
  gamma_est=mean(result_dep$theta)
  
  
  main_predict2=main_predict2
  #main_predict2$x1=main_predict2$mi_st
  #main_predict2$x2=main_predict2$GEND01.x
  #main_predict2$x3=main_predict2$RACE01.x
  #main_predict2$x4=main_predict2$AVZMSYS
  
  
  
  n_patient<-length(unique(main_predict2$id)) 
  x1<-main_predict2$x1
  x2<-main_predict2$x2
  x3<-main_predict2$x3
  x4<-main_predict2$x4
  
 
  
  X<-unique(main_predict2$id)
  for (i in 1:length(X)){
    temp<-main_predict2[main_predict2$id==X[i],]
    temp1<-temp$y1
    main_predict2[main_predict2$id==X[i],"maxni"]<-rep(length(temp1),length(temp1))
  }
  J<- max(main_predict2$maxni)
  y1=main_predict2$y1
  s1=main_predict2$s1
  ########################################################################
  ##########begin prediction and plot the curve of the Brier Scores ########
  ##########################################################################
  ##########prediction using SFM #################
  
  
  
  n_predict<-length(unique(main_predict2$id))
  
  t_s=c(0.06,0.09,0.12)
  bs_t<-list()
  bs_t[[1]]<-matrix(c(rep(0,30)),nrow=15,ncol=2)
  bs_t[[2]]<-matrix(c(rep(0,24)),nrow=12,ncol=2)
  bs_t[[3]]<-matrix(c(rep(0,18)),nrow=9,ncol=2)
  bs_t_terminal<-list()
  bs_t_terminal[[1]]<-matrix(c(rep(0,30)),nrow=15,ncol=2)
  bs_t_terminal[[2]]<-matrix(c(rep(0,24)),nrow=12,ncol=2)
  bs_t_terminal[[3]]<-matrix(c(rep(0,18)),nrow=9,ncol=2)
  
  for(t_prime in t_s){
    
    t_f=seq(t_prime,0.2,0.01)
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
        if(nrow(main_predict_i2)==1| main_predict_i2$y1[1]<t_prime|is.na(main_predict_i2$y1[1])) {next}
        main_predict_p2=main_predict_i2[main_predict_i2$t3<t_prime,]
        if(length(main_predict_p2$y1)==0) {next}
        s1_t2[k]= s1[k]
        y1_t2[k]= y1[k]
        s1_t[k]= s1[k]
        y1_t[k]= y1[k]
        main_predict_p2$y1=t_prime
        main_predict_p2$y2[nrow(main_predict_p2)]=t_prime-main_predict_p2$t3[nrow(main_predict_p2)]
        main_predict_p2$s2[nrow(main_predict_p2)]=0
        main_predict_s2=main_predict_i2[nrow(main_predict_i2),]
        
        set.seed(n.sim*100) 
        
        
        b1<-matrix(cbind(result_dep$b1,result_dep$b2,result_dep$b3,result_dep$b4),nrow=n_sample,ncol=4,byrow = T)
        b2<-matrix(cbind(result_dep$b5,result_dep$b6,result_dep$b7,result_dep$b8),nrow=n_sample,ncol=4,byrow = T)
        x<-matrix(cbind(main_predict_p2$x1,main_predict_p2$x2,main_predict_p2$x3,main_predict_p2$x4),nrow=length(main_predict_p2$x1),ncol=4,byrow = T)
        
        
        phiN_sample_post=phiN_sample=c()
        phiN_sample[1]=0
        phiN_sample_post[1]=0
        set.seed(n.sim*100)
        j=1
        for(m in slice){
          print(m)
          #####MCMC begin######
          for(un_i in 2:10){
            #matrix<-as.matrix(result$theta,nrow=1000,ncol=length(result$theta)/1000,byrow=F)
            phiN_sample[un_i]<- metroplis_phi_id_dep(b1[m,], b2[m,], result_dep$b01[m], result_dep$b02[m],
                                                     result_dep$sigma[m],result_dep$mu[m],result_dep$theta[m],x, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s1,main_predict_p2$s2,
                                                     nrow(main_predict_p2),phiN_sample[un_i-1],0.1,p=dim(x)[2])
            
          }
          phiN_sample_post[j]=phiN_sample[10] ###keep the final random effect sample into the random effects for thet(m)###
          j=j+1  
        }
        
        
        phiN_est2=mean(phiN_sample_post)
        predict2[k]=exp(-t_final*exp(b1_est*main_predict_p2$x1[1]+b2_est*main_predict_p2$x2[1]+b3_est*main_predict_p2$x3[1]+b4_est*main_predict_p2$x4[1]+b01_est+phiN_est2))  
        #predict_terminal[k]=exp(-t_final*exp(b1_est*main_predict_p$x1[1]+b2_est*main_predict_p$x2[1]+b01_est+phiN_est))  
        if (is.na(y1_t2[k])) {y1_t2[k]<-rexp(1,rate=1/(b1_est*main_predict_p2$x1[1]+b2_est*main_predict_p2$x2[1]+b3_est*main_predict_p2$x3[1]+b4_est*main_predict_p2$x4[1]+b01_est+phiN_est2))} 
        if (is.na(s1_t2[k])) {s1_t2[k]<-ifelse(y1_t2<=t_final,1,0)} 
        
        
        b1<-matrix(cbind(result_cop$b1,result_cop$b2,result_cop$b3,result_cop$b4),nrow=n_sample,ncol=4,byrow = T)
        b2<-matrix(cbind(result_cop$b5,result_cop$b6,result_cop$b7,result_cop$b8),nrow=n_sample,ncol=4,byrow = T)
        phiN_sample_post=phiN_sample=c()
        phiN_sample[1]=0
        phiN_sample_post[1]=0
        set.seed(n.sim*100)
        j=1
        for(m in slice){
          print(m)
          #####MCMC begin######
          for(un_i in 2:10){
            
            phiN_sample[un_i]<- metroplis_phi_id(b1[m,], b2[m,], result_cop$b01[m], result_cop$b02[m],
                                                 result_cop$sigma[m],result_cop$theta[m],x, main_predict_p2$y1, main_predict_p2$y2, main_predict_p2$s1,main_predict_p2$s2,
                                                 nrow(main_predict_p2),phiN_sample[un_i-1],0.1,p=dim(x)[2])
            
          }
          phiN_sample_post[j]=phiN_sample[10] ###keep the final random effect sample into the random effects for thet(m)###
          j=j+1  
        }
        
        phiN_est=mean(phiN_sample_post)
        
        predict_terminal[k]=exp(-t_final*exp(b1_est1*main_predict_p2$x1[1]+b2_est1*main_predict_p2$x2[1]+b3_est1*main_predict_p2$x3[1]+b4_est1*main_predict_p2$x4[1]+b01_est1+phiN_est))  
        
        
        if (is.na(y1_t[k])) {y1_t[k]<-rexp(1,rate=1/(b1_est1*main_predict_p2$x1[1]+b2_est1*main_predict_p2$x2[1]+b3_est1*main_predict_p2$x3[1]+b4_est1*main_predict_p2$x4[1]+b01_est1+phiN_est))} 
        if (is.na(s1_t[k])) {s1_t[k]<-ifelse(y1_t<=t_final,1,0)} 
        
        }
      
      
      
      predict2=na.omit(predict2)
      
      censor2=1-na.omit(s1_t2)
      s1_t2=na.omit(s1_t2)
      y1_t2=na.omit(y1_t2)
      
      predict_terminal=na.omit(predict_terminal)
      censor=1-na.omit(s1_t)
      s1_t=na.omit(s1_t)
      y1_t=na.omit(y1_t)
      #predict_time=seq(15000,max(y1_t2),1000)
      #predict_time2=seq(0.1,max(y1_t2),0.01)
      #s_score_terminal=sbrier(Surv(y1_t,s1_t),as.vector(predict_terminal[1:length(y1_t)]),btime=t_final)
      #s_score_terminal=sbrier(Surv(y1_t,s1_t),as.vector(predict_terminal)[1:length(y1_t)],btime=t_final)
      s_score_terminal=sbrier(Surv(y1_t,s1_t),as.vector(predict_terminal),btime=t_final)
      cat(capture.output(print(predict_terminal), file=paste0("predict_terminal",t_final,".txt")))
      #s_score=sbrier(Surv(y1_t2,s1_t2),as.vector(predict2[1:length(y1_t2)]),btime =t_final)
      #s_score=sbrier(Surv(y1_t2,s1_t2),as.vector(predict2)[1:length(y1_t2)],btime =t_final)
      s_score=sbrier(Surv(y1_t2,s1_t2),as.vector(predict2),btime =t_final)
      cat(capture.output(print(predict2), file=paste0("predict2",t_final,".txt")))
      bs_t_terminal[[idd]]=rbind(c(t_final,s_score_terminal),bs_t_terminal[[idd]])
      bs_t[[idd]]=rbind(c(t_final,s_score),bs_t[[idd]])
    }
    
  }
  
  return(list(bs_t,bs_t_terminal))
  
}



length(unique(Stroke$id)) 


######## croaa-validation #####
numCores <- 4
registerDoParallel(numCores)
nsim<-40
start_time <- Sys.time()
result<- foreach(n.sim=1:nsim,.errorhandling = "remove") %dopar% {
  set.seed(n.sim*100)
  ind<-sample(1:length(unique(Stroke$id)),1000,replace = FALSE)
  work.data<-Stroke[Stroke$id %in% ind,]
  train_ind<-sample(1:length(unique(work.data$id)),600,replace = FALSE)
  train.data<-work.data[work.data$id %in% train_ind,]
  test.data<-work.data[work.data$id %in% unique(work.data$id)[-train_ind],]
  test.data<-test.data[order(test.data$id),]
  X<-unique(test.data$id)
  id_new<-list()
  for (i in 1 :length(X)){
    tmp<-test.data[test.data$id==X[i],]
    id_new[[i]]<-rep(i,nrow(tmp))
  }
  
  test.data$id_new<-unlist(id_new)
  test.data$id<-test.data$id_new
  train.data<-train.data[order(train.data$id),]
  X<-unique(train.data$id)
  id_new<-list()
  for (i in 1 :length(X)){
    tmp<-train.data[train.data$id==X[i],]
    id_new[[i]]<-rep(i,nrow(tmp))
  }
  
  train.data$id_new<-unlist(id_new)
  train.data$id<-train.data$id_new
  
  test.data<-test.data[order(test.data$id),]
  train.data<-train.data[order(train.data$id),]
  train.data$y1<-train.data$y1/(365.25*100)
  train.data$y2<-train.data$y2/(365.25*100)
  train.data$t3<-train.data$t3/(365.25*100)
  test.data$y1<-test.data$y1/(365.25*100)
  test.data$y2<-test.data$y2/(365.25*100)
  test.data$t3<-test.data$t3/(365.25*100)
  #train.data$x1<-scale(train.data$x1,scale = T)
  #train.data$x1<-as.numeric(train.data$x1)
  #test.data$x1<-scale(test.data$x1,scale = T)
  #test.data$x1<-as.numeric(test.data$x1)
  ##### generate initial values for parameters ####
  
  ### call function to return brier score ####
  
  out(train.data,test.data,n.sim = n.sim)
  
}

end_time <- Sys.time()
end_time - start_time
setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/02data")
saveRDS(result,paste0("CHS_st_",nsim,"_0119.rds"))
result<- readRDS(file.choose())
############# function for plot ###########
BS<-function(result){
  
  
  
  bs_t=list()
  bs_t_terminal=list()
  
  bs_t_up=list()
  bs_t_terminal_up=list()
  
  bs_t_sd=list()
  bs_t_terminal_sd=list()
  
  bs_t[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_sd[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t[[1]][,2]<-c(rep(0,10))
  bs_t[[1]][,1]<-c(rep(0,10))
  bs_t[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_sd[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t[[2]][,2]<-c(rep(0,7))
  bs_t[[2]][,1]<-c(rep(0,7))
  bs_t[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_sd[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t[[3]][,2]<-c(rep(0,4))
  bs_t[[3]][,1]<-c(rep(0,4))
  bs_t_terminal[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal_sd[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal[[1]][,2]<-c(rep(0,10))
  bs_t_terminal[[1]][,1]<-c(rep(0,10))
  bs_t_terminal[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal_sd[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal[[2]][,2]<-c(rep(0,7))
  bs_t_terminal[[2]][,1]<-c(rep(0,7))
  bs_t_terminal[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal_sd[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal[[3]][,2]<-c(rep(0,4))
  bs_t_terminal[[3]][,1]<-c(rep(0,4))
  
  
  
  bs_t_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  bs_t_terminal_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_terminal_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_terminal_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  
  
  for (i in 1:length(result)){
    for (j in 1:10) {
      temp1<-ifelse(is.null(result[[i]][[1]][[1]][j,2]),0,result[[i]][[1]][[1]][j,2])
      bs_t_up[[1]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[1]][j,2]),0,result[[i]][[2]][[1]][j,2])
      bs_t_terminal_up[[1]][i,j]<-temp2
    }
  }
  
  
  for (i in 1:length(result)){
    for (j in 1:7) {
      temp1<-ifelse(is.null(result[[i]][[1]][[2]][j,2]),0,result[[i]][[1]][[2]][j,2])
      bs_t_up[[2]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[2]][j,2]),0,result[[i]][[2]][[2]][j,2])
      bs_t_terminal_up[[2]][i,j]<-temp2
    }
  }
  for (i in 1:length(result)){
    for (j in 1:4) {
      temp1<-ifelse(is.null(result[[i]][[1]][[3]][j,2]),0,result[[i]][[1]][[3]][j,2])
      bs_t_up[[3]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[3]][j,2]),0,result[[i]][[2]][[3]][j,2])
      bs_t_terminal_up[[3]][i,j]<-temp2
    }
  }
  
  #bs_t_sd[[1]][1:15,2]<-apply(bs_t_up[[1]], 2, sd)
  # bs_t_terminal_sd[[1]][1:15,2]<-apply(bs_t_terminal_up[[1]], 2, sd)
  #bs_t_sd[[2]][1:9,2]<-apply(bs_t_up[[2]], 2, sd)
  #bs_t_terminal_sd[[2]][1:5,2]<-apply(bs_t_terminal_up[[2]], 2, sd)
  #  bs_t_sd[[3]][1:3,2]<-apply(bs_t_up[[3]], 2, sd)
  #bs_t_terminal_sd[[3]][1:2,2]<-apply(bs_t_terminal_up[[3]], 2, sd)
  
  bs_t[[1]][1:10,2]<-apply(bs_t_up[[1]], 2, mean)
  bs_t_terminal[[1]][1:10,2]<-apply(bs_t_terminal_up[[1]], 2, mean)
  bs_t[[2]][1:7,2]<-apply(bs_t_up[[2]], 2, mean)
  bs_t_terminal[[2]][1:7,2]<-apply(bs_t_terminal_up[[2]], 2, mean)
  bs_t[[3]][1:4,2]<-apply(bs_t_up[[3]], 2, mean)
  bs_t_terminal[[3]][1:4,2]<-apply(bs_t_terminal_up[[3]], 2, mean)
  
  
  bs_t[[1]][1:10,1]<-seq(0.06,0.15,0.01)
  bs_t[[2]][1:7,1]<-seq(0.09,0.15,0.01)
  bs_t[[3]][1:4,1]<-seq(0.12,0.15,0.01)
  #bs_t_sd[[1]][1:15,1]<-c(3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000)
  #bs_t_sd[[2]][1:9,1]<-c(6000,6500,7000,7500,8000,8500,9000,9500,10000)
  #bs_t_sd[[3]][1:3,1]<-c(9000,9500,10000)
  
  
  bs_t_terminal[[1]][1:10]<-seq(0.06,0.15,0.01)
  bs_t_terminal[[2]][1:7]<-seq(0.09,0.15,0.01)
  bs_t_terminal[[3]][1:4]<-seq(0.12,0.15,0.01)
  #bs_t_terminal_sd[[1]][1:8]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  #bs_t_terminal_sd[[2]][1:5]<-c(0.1,0.09,0.08,0.07,0.06)
  #bs_t_terminal_sd[[3]][1:2]<-c(0.1,0.09)
  
  setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03results/")
  #setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03preresults2/")
  #res<-list(bs_t,bs_t_terminal,bs_t_sd,bs_t_terminal_sd)
  res<-list(bs_t,bs_t_sd)
  #sink(paste0("output",idx,".txt"))
  #print(res)
  #sink()
  cat(capture.output(print(res), file="output_CHS_Stroke_0119.txt"))
  #cat(capture.output(print(bs_t), file=paste0("BS_terminal",idx,".txt")))
  #cat(capture.output(print(bs_t_terminal), file=paste0("BS_terminal_old",idx,".txt")))
  #cat(capture.output(print(bs_t_sd), file=paste0("BS_terminal_SD",idx,".txt")))
  #cat(capture.output(print(bs_t_terminal_sd), file=paste0("BS_terminal_SD_old",idx,".txt")))
  
  #temp<-as.data.frame(cbind(bs_t[[1]],bs_t_sd[[1]][,2],bs_t_terminal[[1]][,2],
  #                         bs_t_terminal_sd[[1]][,2]))
  #colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  #write.csv(temp, file=paste0("BS_terminal_0.03_",idx,"_.csv"))
  temp<-as.data.frame(cbind(bs_t[[1]],bs_t_terminal[[1]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_mean_old")
  write.csv(temp, file="BS_terminal_6year_CHS_Stroke_0119.csv")
  
  # temp<-as.data.frame(cbind(bs_t[[2]],bs_t_sd[[2]][,2],bs_t_terminal[[2]][,2],
  #                           bs_t_terminal_sd[[2]][,2]))
  # colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  # write.csv(temp, file=paste0("BS_terminal_0.06_",idx,"_.csv"))
  
  temp<-as.data.frame(cbind(bs_t[[2]],bs_t_terminal[[2]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_mean_old")
  write.csv(temp, file="BS_terminal_9year_CHS_Stroke_0119.csv")
  
  
  #temp<-as.data.frame(cbind(bs_t[[3]],bs_t_sd[[3]][,2],bs_t_terminal[[3]][,2],
  #                          bs_t_terminal_sd[[3]][,2]))
  #colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  #write.csv(temp, file=paste0("BS_terminal_0.09_",idx,"_.csv"))
  
  temp<-as.data.frame(cbind(bs_t[[3]],bs_t_terminal[[3]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_mean_old")
  write.csv(temp, file="BS_terminal_12year_CHS_Stroke_0119.csv")
  
  
  png("BS_CHS_Stroke_0119.png",width=900, height=900)
  
  t_s=c(6,9,12)
  
  for(t_prime in t_s){
    
    idxx=which(t_s==t_prime)
    
    
    if(t_prime==t_s[1]) {
      
      bs_t_terminal[[1]][,1]=bs_t_terminal[[1]][,1]*100
      bs_t_terminal[[2]][,1]=bs_t_terminal[[2]][,1]*100
      bs_t_terminal[[3]][,1]=bs_t_terminal[[3]][,1]*100
      
      bs_t[[1]][,1]=bs_t[[1]][,1]*100
      bs_t[[2]][,1]=bs_t[[2]][,1]*100
      bs_t[[3]][,1]=bs_t[[3]][,1]*100
      
      plot(bs_t[[1]][1:10],bs_t[[1]][,2],xlab="",type="o",xlim=c(5,15),ylim=c(0,1),main="",lty=1,col = "red", ylab = "",pch=16)
      #arrows(x0=bs_t[[1]][1:8],y0=bs_t[[1]][,2]-bs_t_sd[[1]][,2], x1=bs_t[[1]][1:8], y1=bs_t[[1]][,2]+bs_t_sd[[1]][,2], code=3, angle=90, length=0.1,col = "red")
      legend("topleft",c("t'=6-yr","t'=9-yr","t'=12-yr"),lty=1:3,lwd=1,pch=1,col = "red",title= "Time-varying JFCM",cex = 2)
      lines(bs_t_terminal[[1]],type="o",lty=1)
      mtext(side=1, line=2, "year", font=2.5,cex=1.2)
      mtext(side=2, line=3, "Brier score",font=2.5,cex=1.2)
      #arrows(x0=bs_t_terminal[[1]][1:8],y0=bs_t_terminal[[1]][,2]-bs_t_terminal_sd[[1]][,2], x1=bs_t_terminal[[1]][1:8], y1=bs_t_terminal[[1]][,2]+bs_t_terminal_sd[[1]][,2], code=3, angle=90, length=0.1)
    }
    else {
      lines(bs_t[[2]],type="o",lty=2,col = "red")
      #arrows(x0=bs_t[[2]][1:5],y0=bs_t[[2]][,2]-bs_t_sd[[2]][,2], x1=bs_t[[2]][1:5], y1=bs_t[[2]][,2]+bs_t_sd[[2]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t[[3]],type="o",lty=3,col = "red")
      #arrows(x0=bs_t[[3]][1:2],y0=bs_t[[3]][,2]-bs_t_sd[[3]][,2], x1=bs_t[[3]][1:2], y1=bs_t[[3]][,2]+bs_t_sd[[3]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t_terminal[[2]],type="o",lty=2)
      #arrows(x0=bs_t_terminal[[2]][1:5],y0=bs_t_terminal[[2]][,2]-bs_t_terminal_sd[[2]][,2], x1=bs_t_terminal[[2]][1:5], y1=bs_t_terminal[[2]][,2]+bs_t_terminal_sd[[2]][,2], code=3, angle=90, length=0.1)
      lines(bs_t_terminal[[3]],type="o",lty=3)
      #arrows(x0=bs_t_terminal[[3]][1:2],y0=bs_t_terminal[[3]][,2]-bs_t_terminal_sd[[3]][,2], x1=bs_t_terminal[[3]][1:2], y1=bs_t_terminal[[3]][,2]+bs_t_terminal_sd[[3]][,2], code=3, angle=90, length=0.1)
    }
    abline(v=t_prime,lty=2)
    if(t_prime==t_s[1]) {axis(side=1,seq(5,15,1))
      axis(side=2, seq(0,1))}
  }
  
  legend("bottomleft",c("t'=6-yr","t'=9-yr","t'=12-yr"),lty=1:3,lwd=1,pch=1,title= "JFCM",cex = 2)
  ###################################################################################
  
  dev.off()
  
  #return(list(bs_t, bs_t_terminal,bs_t_sd,bs_t_terminal_sd))
  ###################################################################################
  
}



BS(result)



out2<-function(st_final,n.sim){
  
  set.seed(n.sim*100)
  ind<-sample(1:length(unique(st_final$id)),1000,replace = FALSE)
  work.data<-st_final[st_final$id %in% ind,]
  train_ind<-sample(1:length(unique(work.data$id)),600,replace = FALSE)
  train.data<-work.data[work.data$id %in% train_ind,]
  test.data<-work.data[work.data$id %in% unique(work.data$id)[-train_ind],]
  test.data<-test.data[order(test.data$id),]
  X<-unique(test.data$id)
  id_new<-list()
  for (i in 1 :length(X)){
    tmp<-test.data[test.data$id==X[i],]
    id_new[[i]]<-rep(i,nrow(tmp))
  }
  
  test.data$id_new<-unlist(id_new)
  test.data$id<-test.data$id_new
  train.data<-train.data[order(train.data$id),]
  X<-unique(train.data$id)
  id_new<-list()
  for (i in 1 :length(X)){
    tmp<-train.data[train.data$id==X[i],]
    id_new[[i]]<-rep(i,nrow(tmp))
  }
  
  train.data$id_new<-unlist(id_new)
  train.data$id<-train.data$id_new
  
  test.data<-test.data[order(test.data$id),]
  train.data<-train.data[order(train.data$id),]
  train.data$y1<-train.data$y1/(365.25*100)
  train.data$y2<-train.data$y2/(365.25*100)
  train.data$t3<-train.data$t3/(365.25*100)
  test.data$y1<-test.data$y1/(365.25*100)
  test.data$y2<-test.data$y2/(365.25*100)
  test.data$t3<-test.data$t3/(365.25*100)
  
  
  modJoint.gap<- frailtyPenal(Surv(y2,s2) ~ cluster(id) +
                                x1+x2+x3+x4+terminal(s1),
                              formula.terminalEvent = ~ x1+x2+x3+x4,data = train.data, n.knots = 8, kappa = c(2,3))
  
  
  b1_est_SF<-modJoint.gap$coef[1]
  b2_est_SF<-modJoint.gap$coef[2]
  b3_est_SF<-modJoint.gap$coef[3]
  b4_est_SF<-modJoint.gap$coef[4]
  modJoint.gap<- frailtyPenal(Surv(y2,s2) ~ cluster(id) +
                                x1+x2+x3+x4+terminal(s1),
                              formula.terminalEvent = ~ x1+x2+x3+x4,data = test.data, n.knots = 8, kappa = c(2,3))
  main_predict2<-test.data
  phi_est_SF<-modJoint.gap$alpha
  y1<-vector()
  s1<-vector()
  prediction<-vector()
  for (i in unique(main_predict2$id)){
    ni<-dim(main_predict2[main_predict2$id==i,])[1]
    y1[i]<-max(main_predict2[main_predict2$id==i,]$y1)
    s1[i]<-main_predict2[main_predict2$id==i,]$s1[ni]
    
  }
  t_s=c(0.06,0.09,0.12)
  bs_t_terminal_SF<-list()
  bs_t_terminal_SF[[1]]<-matrix(c(rep(0,30)),nrow=15,ncol=2)
  bs_t_terminal_SF[[2]]<-matrix(c(rep(0,24)),nrow=12,ncol=2)
  bs_t_terminal_SF[[3]]<-matrix(c(rep(0,18)),nrow=9,ncol=2)
  for(t_prime in t_s){
    
    t_f=seq(t_prime,0.2,0.01)
    idd<-which(t_s==t_prime)
    for (t_final in t_f){
      phi_est_SF=rnorm(1,phi_est_SF,0.1)
      for(k in 1:length(unique(main_predict2$id))){
        
        main_predict_p2=main_predict2[main_predict2$id==k,]
        
        ############## using existing function in frailtypack package ####################
        set.seed(n.sim*100)
        
        
        # pred<-prediction(modJoint.gap, main_predict, t =t_prime-0.01, window =seq(0.01, 0.1, 0.01), MC.sample = 60)
        prediction[k]<-exp(-t_final*exp(b1_est_SF*main_predict_p2$x1[1]+b2_est_SF*main_predict_p2$x2[1]+
                                          b3_est_SF*main_predict_p2$x3[1]+b4_est_SF*main_predict_p2$x4[1]+phi_est_SF))  
      }
      
      
      
      s_score_terminal_SF=sbrier(Surv(y1,s1),prediction,btime =t_final)
      bs_t_terminal_SF[[idd]]=rbind(c(t_final,s_score_terminal_SF),bs_t_terminal_SF[[idd]])
      
    }
  }
  
  
  return(list(bs_t_terminal_SF))
}

res<-foreach(n.sim=1:nsim,.errorhandling = "remove") %dopar% {out2(Stroke,n.sim=n.sim)}


BS<-function(result,res){
  
  
  
  
  bs_t=list()
  bs_t_terminal=list()
  bs_t_terminal_SF=list()
  
  bs_t_up=list()
  bs_t_terminal_up=list()
  bs_t_terminal_SF_up=list()
  bs_t_sd=list()
  bs_t_terminal_sd=list()
  bs_t_terminal_SF_sd=list()
  
  
  
  
  
  
  bs_t[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_sd[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t[[1]][,2]<-c(rep(0,10))
  bs_t[[1]][,1]<-c(rep(0,10))
  bs_t[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_sd[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t[[2]][,2]<-c(rep(0,7))
  bs_t[[2]][,1]<-c(rep(0,7))
  bs_t[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_sd[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t[[3]][,2]<-c(rep(0,4))
  bs_t[[3]][,1]<-c(rep(0,4))
  bs_t_terminal[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal_sd[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal[[1]][,2]<-c(rep(0,10))
  bs_t_terminal[[1]][,1]<-c(rep(0,10))
  bs_t_terminal[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal_sd[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal[[2]][,2]<-c(rep(0,7))
  bs_t_terminal[[2]][,1]<-c(rep(0,7))
  bs_t_terminal[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal_sd[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal[[3]][,2]<-c(rep(0,4))
  bs_t_terminal[[3]][,1]<-c(rep(0,4))
  
  
  
  bs_t_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  bs_t_terminal_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_terminal_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_terminal_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  
  
  for (i in 1:length(result)){
    for (j in 1:10) {
      temp1<-ifelse(is.null(result[[i]][[1]][[1]][j,2]),0,result[[i]][[1]][[1]][j,2])
      bs_t_up[[1]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[1]][j,2]),0,result[[i]][[2]][[1]][j,2])
      bs_t_terminal_up[[1]][i,j]<-temp2
    }
  }
  
  
  for (i in 1:length(result)){
    for (j in 1:7) {
      temp1<-ifelse(is.null(result[[i]][[1]][[2]][j,2]),0,result[[i]][[1]][[2]][j,2])
      bs_t_up[[2]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[2]][j,2]),0,result[[i]][[2]][[2]][j,2])
      bs_t_terminal_up[[2]][i,j]<-temp2
    }
  }
  for (i in 1:length(result)){
    for (j in 1:4) {
      temp1<-ifelse(is.null(result[[i]][[1]][[3]][j,2]),0,result[[i]][[1]][[3]][j,2])
      bs_t_up[[3]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[3]][j,2]),0,result[[i]][[2]][[3]][j,2])
      bs_t_terminal_up[[3]][i,j]<-temp2
    }
  }
  
  bs_t_terminal_SF[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal_SF_sd[[1]]<-matrix(c(rep(0,20)),nrow=10,ncol=2)
  bs_t_terminal_SF[[1]][,2]<-c(rep(0,10))
  bs_t_terminal_SF[[1]][,1]<-c(rep(0,10))
  bs_t_terminal_SF[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal_SF_sd[[2]]<-matrix(c(rep(0,14)),nrow=7,ncol=2)
  bs_t_terminal_SF[[2]][,2]<-c(rep(0,7))
  bs_t_terminal_SF[[2]][,1]<-c(rep(0,7))
  bs_t_terminal_SF[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal_SF_sd[[3]]<-matrix(c(rep(0,8)),nrow=4,ncol=2)
  bs_t_terminal_SF[[3]][,2]<-c(rep(0,4))
  bs_t_terminal_SF[[3]][,1]<-c(rep(0,4))
  
  
  
  for (i in 1:length(res)) {
    bs_t_terminal_SF[[1]]<-res[[i]][[1]][[1]][1:10,]+bs_t_terminal_SF[[1]][1:10]
    bs_t_terminal_SF[[2]]<-res[[i]][[1]][[2]][1:7,]+bs_t_terminal_SF[[2]][1:7]
    bs_t_terminal_SF[[3]]<-res[[i]][[1]][[3]][1:4,]+bs_t_terminal_SF[[3]][1:4]
  }
  
  bs_t_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  bs_t_terminal_up[[1]]<-matrix(c(rep(0,10*length(result))),nrow=length(result),ncol=10)
  bs_t_terminal_up[[2]]<-matrix(c(rep(0,7*length(result))),nrow=length(result),ncol=7)
  bs_t_terminal_up[[3]]<-matrix(c(rep(0,4*length(result))),nrow=length(result),ncol=4)
  
  bs_t_terminal_SF_up[[1]]<-matrix(c(rep(0,10*length(res))),nrow=length(res),ncol=10)
  bs_t_terminal_SF_up[[2]]<-matrix(c(rep(0,7*length(res))),nrow=length(res),ncol=7)
  bs_t_terminal_SF_up[[3]]<-matrix(c(rep(0,4*length(res))),nrow=length(res),ncol=4)
  
  
  for (i in 1:length(result)){
    for (j in 1:10) {
      temp1<-ifelse(is.null(result[[i]][[1]][[1]][j,2]),0,result[[i]][[1]][[1]][j,2])
      bs_t_up[[1]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[1]][j,2]),0,result[[i]][[2]][[1]][j,2])
      bs_t_terminal_up[[1]][i,j]<-temp2
    }
  }
  
  
  for (i in 1:length(result)){
    for (j in 1:7) {
      temp1<-ifelse(is.null(result[[i]][[1]][[2]][j,2]),0,result[[i]][[1]][[2]][j,2])
      bs_t_up[[2]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[2]][j,2]),0,result[[i]][[2]][[2]][j,2])
      bs_t_terminal_up[[2]][i,j]<-temp2
    }
  }
  for (i in 1:length(result)){
    for (j in 1:4) {
      temp1<-ifelse(is.null(result[[i]][[1]][[3]][j,2]),0,result[[i]][[1]][[3]][j,2])
      bs_t_up[[3]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[3]][j,2]),0,result[[i]][[2]][[3]][j,2])
      bs_t_terminal_up[[3]][i,j]<-temp2
    }
  }
  
  bs_t_sd[[1]][1:10,2]<-apply(bs_t_up[[1]], 2, sd)
  bs_t_terminal_sd[[1]][1:10,2]<-apply(bs_t_terminal_up[[1]], 2, sd)
  bs_t_sd[[2]][1:7,2]<-apply(bs_t_up[[2]], 2, sd)
  bs_t_terminal_sd[[2]][1:7,2]<-apply(bs_t_terminal_up[[2]], 2, sd)
  bs_t_sd[[3]][1:4,2]<-apply(bs_t_up[[3]], 2, sd)
  bs_t_terminal_sd[[3]][1:4,2]<-apply(bs_t_terminal_up[[3]], 2, sd)
  
  bs_t[[1]][1:10,2]<-apply(bs_t_up[[1]], 2, mean)
  bs_t_terminal[[1]][1:10,2]<-apply(bs_t_terminal_up[[1]], 2, mean)
  bs_t[[2]][1:7,2]<-apply(bs_t_up[[2]], 2, mean)
  bs_t_terminal[[2]][1:7,2]<-apply(bs_t_terminal_up[[2]], 2, mean)
  bs_t[[3]][1:4,2]<-apply(bs_t_up[[3]], 2, mean)
  bs_t_terminal[[3]][1:4,2]<-apply(bs_t_terminal_up[[3]], 2, mean)
  
  
  bs_t[[1]][1:10,1]<-seq(0.06,0.15,0.01)
  bs_t[[2]][1:7,1]<-seq(0.09,0.15,0.01)
  bs_t[[3]][1:4,1]<-seq(0.12,0.15,0.01)
  
  
  bs_t_terminal[[1]][1:10]<-seq(0.06,0.15,0.01)
  bs_t_terminal[[2]][1:7]<-seq(0.09,0.15,0.01)
  bs_t_terminal[[3]][1:4]<-seq(0.12,0.15,0.01)
  
  for (i in 1:length(res)){
    for (j in 1:10) {
      bs_t_terminal_SF_up[[1]][i,j]<-res[[i]][[1]][[1]][j,2]
    }
  }
  
  for (i in 1:length(res)){
    for (j in 1:7) {
      bs_t_terminal_SF_up[[2]][i,j]<-res[[i]][[1]][[2]][j,2]
    }
  }
  for (i in 1:length(res)){
    for (j in 1:4) {
      bs_t_terminal_SF_up[[3]][i,j]<-res[[i]][[1]][[3]][j,2]
    }
  }
  bs_t_terminal_SF_up[[1]]<-na.omit(bs_t_terminal_SF_up[[1]])
  bs_t_terminal_SF_up[[2]]<-na.omit(bs_t_terminal_SF_up[[2]])
  bs_t_terminal_SF_up[[3]]<-na.omit(bs_t_terminal_SF_up[[3]])
  
  bs_t_sd[[1]][1:10,2]<-apply(bs_t_up[[1]], 2, sd)
  bs_t_terminal_sd[[1]][1:10,2]<-apply(bs_t_terminal_up[[1]], 2, sd)
  bs_t_terminal_SF_sd[[1]][1:10,2]<-apply(bs_t_terminal_SF_up[[1]], 2, sd)
  bs_t_terminal_SF[[1]][1:10,2]<-apply(bs_t_terminal_SF_up[[1]], 2, mean)
  bs_t_sd[[2]][1:7,2]<-apply(bs_t_up[[2]], 2, sd)
  bs_t_terminal_sd[[2]][1:7,2]<-apply(bs_t_terminal_up[[2]], 2, sd)
  bs_t_terminal_SF_sd[[2]][1:7,2]<-apply(bs_t_terminal_SF_up[[2]], 2, sd)
  bs_t_terminal_SF[[2]][1:7,2]<-apply(bs_t_terminal_SF_up[[2]], 2, mean)
  bs_t_sd[[3]][1:4,2]<-apply(bs_t_up[[3]], 2, sd)
  bs_t_terminal_sd[[3]][1:4,2]<-apply(bs_t_terminal_up[[3]], 2, sd)
  bs_t_terminal_SF_sd[[3]][1:4,2]<-apply(bs_t_terminal_SF_up[[3]], 2, sd)
  bs_t_terminal_SF[[3]][1:4,2]<-apply(bs_t_terminal_SF_up[[3]], 2, mean)
  
  #bs_t_terminal_SF[[1]]<-bs_t_terminal_SF[[1]]/length(res)
  #bs_t_terminal_SF[[2]]<-bs_t_terminal_SF[[2]]/length(res)
  #bs_t_terminal_SF[[3]]<-bs_t_terminal_SF[[3]]/length(res)
  bs_t_terminal_SF[[1]][1:10]<-seq(0.06,0.15,0.01)
  bs_t_terminal_SF[[2]][1:7]<-seq(0.09,0.15,0.01)
  bs_t_terminal_SF[[3]][1:4]<-seq(0.12,0.15,0.01)
  bs_t_terminal_SF_sd[[1]][1:10]<-seq(0.06,0.15,0.01)
  bs_t_terminal_SF_sd[[2]][1:7]<-seq(0.09,0.15,0.01)
  bs_t_terminal_SF_sd[[3]][1:4]<-seq(0.12,0.15,0.01)
  
  setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03results/")
  #setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03preresults2/")
  res<-list(bs_t,bs_t_terminal,bs_t_terminal_SF,bs_t_sd,bs_t_terminal_sd,bs_t_terminal_SF_sd)
  #sink(paste0("output",idx,".txt"))
  #print(res)
  #sink()
  cat(capture.output(print(res), file=paste0("output_WSFM_CHS_Stroke_0119.txt")))
  
  temp<-as.data.frame(cbind(bs_t[[1]],bs_t_sd[[1]][,2],bs_t_terminal[[1]][,2],
                            bs_t_terminal_sd[[1]][,2],bs_t_terminal_SF[[1]][,2],
                            bs_t_terminal_SF_sd[[1]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_6yr_WSFM_CHS_Stroke_0119.csv"))
  
  temp<-as.data.frame(cbind(bs_t[[2]],bs_t_sd[[2]][,2],bs_t_terminal[[2]][,2],
                            bs_t_terminal_sd[[2]][,2],bs_t_terminal_SF[[2]][,2],
                            bs_t_terminal_SF_sd[[2]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_9yr_WSFM_CHS_Stroke_0119.csv"))
  
  temp<-as.data.frame(cbind(bs_t[[3]],bs_t_sd[[3]][,2],bs_t_terminal[[3]][,2],
                            bs_t_terminal_sd[[3]][,2],bs_t_terminal_SF[[3]][,2],
                            bs_t_terminal_SF_sd[[3]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_12yr_WSFM_CHS_Stroke_0119.csv"))
  
  
  png(paste0("BS_WSFM_Stroke_CHS_0119.png"),width=900, height=900)
  
  t_s=c(6,9,12)
  
  for(t_prime in t_s){
    
    
    
    idxx=which(t_s==t_prime)
    
    
    if(t_prime==t_s[1]) {
      
      
      bs_t_terminal[[1]][,1]=bs_t_terminal[[1]][,1]*100
      bs_t_terminal[[2]][,1]=bs_t_terminal[[2]][,1]*100
      bs_t_terminal[[3]][,1]=bs_t_terminal[[3]][,1]*100
      
      bs_t_terminal_SF[[1]][,1]=bs_t_terminal_SF[[1]][,1]*100
      bs_t_terminal_SF[[2]][,1]=bs_t_terminal_SF[[2]][,1]*100
      bs_t_terminal_SF[[3]][,1]=bs_t_terminal_SF[[3]][,1]*100
      
      bs_t[[1]][,1]=bs_t[[1]][,1]*100
      bs_t[[2]][,1]=bs_t[[2]][,1]*100
      bs_t[[3]][,1]=bs_t[[3]][,1]*100
      
      plot(bs_t[[1]][1:10],bs_t[[1]][,2],xlab="",type="o",xlim=c(5,15),ylim=c(0,1),main="",lty=1,col = "red", ylab = "",pch=16)
      #arrows(x0=bs_t[[1]][1:8],y0=bs_t[[1]][,2]-bs_t_sd[[1]][,2], x1=bs_t[[1]][1:8], y1=bs_t[[1]][,2]+bs_t_sd[[1]][,2], code=3, angle=90, length=0.1,col = "red")
      legend("topleft",c("t'=6-yr","t'=9-yr","t'=12-yr"),lty=1:3,lwd=1,pch=1,col = "red",title= "Time-varying JFCM",cex = 2)
      lines(bs_t_terminal[[1]],type="o",lty=1)
      lines(bs_t_terminal_SF[[1]],type="o",lty=1,col="blue")
      mtext(side=1, line=2, "year", font=2.5,cex=1.2)
      mtext(side=2, line=3, "Brier score",font=2.5,cex=1.2)
      #arrows(x0=bs_t_terminal[[1]][1:8],y0=bs_t_terminal[[1]][,2]-bs_t_terminal_sd[[1]][,2], x1=bs_t_terminal[[1]][1:8], y1=bs_t_terminal[[1]][,2]+bs_t_terminal_sd[[1]][,2], code=3, angle=90, length=0.1)
    }
    else {
      lines(bs_t[[2]],type="o",lty=2,col = "red")
      #arrows(x0=bs_t[[2]][1:5],y0=bs_t[[2]][,2]-bs_t_sd[[2]][,2], x1=bs_t[[2]][1:5], y1=bs_t[[2]][,2]+bs_t_sd[[2]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t[[3]],type="o",lty=3,col = "red")
      #arrows(x0=bs_t[[3]][1:2],y0=bs_t[[3]][,2]-bs_t_sd[[3]][,2], x1=bs_t[[3]][1:2], y1=bs_t[[3]][,2]+bs_t_sd[[3]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t_terminal[[2]],type="o",lty=2)
      #arrows(x0=bs_t_terminal[[2]][1:5],y0=bs_t_terminal[[2]][,2]-bs_t_terminal_sd[[2]][,2], x1=bs_t_terminal[[2]][1:5], y1=bs_t_terminal[[2]][,2]+bs_t_terminal_sd[[2]][,2], code=3, angle=90, length=0.1)
      lines(bs_t_terminal[[3]],type="o",lty=3)
      #arrows(x0=bs_t_terminal[[3]][1:2],y0=bs_t_terminal[[3]][,2]-bs_t_terminal_sd[[3]][,2], x1=bs_t_terminal[[3]][1:2], y1=bs_t_terminal[[3]][,2]+bs_t_terminal_sd[[3]][,2], code=3, angle=90, length=0.1)
      lines(bs_t_terminal_SF[[2]],type="o",lty=2,col="blue")
      lines(bs_t_terminal_SF[[3]],type="o",lty=3,col="blue")
    }
    abline(v=t_prime,lty=2)
    if(t_prime==t_s[1]) {axis(side=1,seq(5,15,1))
      axis(side=2, seq(0,1))}
  }
  
  legend("top",c("t'=6-yr","t'=9-yr","t'=12-yr"),lty=1:3,lwd=1,pch=1,title= "JFCM",cex = 2)
  legend("topright",c("t'=6-yr","t'=9-yr","t'=12-yr"),lty=1:3,lwd=1,pch=1,title= "SFM",cex = 2,col="blue")
  ###################################################################################
  
  dev.off()
  
  return(list(bs_t, bs_t_terminal,bs_t_terminal_SF,bs_t_sd,bs_t_terminal_sd,bs_t_terminal_SF_sd))
  ###################################################################################
  
}


BS(result,res)
