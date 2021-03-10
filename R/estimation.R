#' ####################################################################
#' #####the function is for log likelihood function of the data########
#' ### function argurments: #######
#' ##b1: the parameter for x1 in the hazard funtion for the terminal event;
#' ##b2: the parameter for x2 in the hazard funtion for the terminal event;
#' ##b3: the parameter for x1 in the hazard funtion for the recurrent event;
#' ##b4: the parameter for x2 in the hazard funtion for the recurrent event;
#' ##b01, b02: baseline hazards for the terminal event and the recurrent events, respectively;
#' ##phi: random effect term denoted as w in paper; {W_i}
#' ##theta: correlation parameter for copula function; {theta_ij}
#' ## Mu: theta_mu average level correlation
#' ##x1, x2: two covariates in the hazard function for the terminal event and the recurrent event;
#' ##y1, y2: observed time for the terminal event and the recurrent event;
#' ##s1, s2: failing indicator for the terminal event and the recurrent event;
#' ##id: subject id;
#' ##N: the sample size of subjects;
#' b01_est1,b02_est1,b1_est1,b2_est1,b3_est1,b4_est1: coefficients in JHCM
#' b01_est,b02_est,b1_est,b2_est,b3_est,b4_est: Coefficiets in TVJHCM
#'####################################################################
#'##########    usage ##############################

#' idx<-4
#'  numCores <- 4
#'  registerDoParallel(numCores)
#'  nsim<-4
#'  result<- foreach(n.sim=1:nsim,.errorhandling = "remove") %dopar% {
#'    Beta_est_sim(idx,n.sim)
#'  }
#'  saveRDS(result,paste0("400_WSF_old_",idx,".rds"))
######################################################
######   generate data and fit the model   ##########
######################################################
#'
#' @param infile Path to the input file
#' @return estimation for coefficient in simulation
#' @export

Beta_est_sim<-function(idx,n.sim){

        n_patient_s=rep(rep(c(100,200,400),each=9),3) ##possible values for sample size, which is set as 100 ,200 or 400##
        theta_s=rep(c(1,1,1,1.5,1.5,1.5,2,2,2),3)              ##possible values for theta, correlation parameter ##########
        theta=theta_s[idx]              ##for scenario idx=4, theta is 1.5##########
        n_patient=n_patient_s[idx]      ##for scenario idx=4, sample size is 100###
        J_s=rep(c(1000,1000,1000),9)    ##for scenario idx=4, J is 100###
        J=J_s[idx]
        UU_s=rep(c(1,0.6,0.4),9)  ## failure rate
        UU=UU_s[idx]
        rho_s=rep(c(0.1,0.3,0.5),each=9)
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
        y1=pmin(t1,c1)    ###observed time####
        s1=(t1<=c1)       ###failing indicator####
        J=J           ###maximum number of the recurrent events###
        uv=rCopula(J,claytonCopula(theta)) ###generate clayton correlated random numbers###
        main=list()


        for(n_i in 1:n_patient){
                t2=0
                repeat{
                        if(s1[n_i]==1){
                                w=runif(1)
                                v=((w^(-theta/(theta+1))-1)*(u[n_i]^(-theta))+1)^(-1/theta)
                                r_new=-log(v)/(lam2[n_i])  ###generate recurrent event time by conditional copula ###
                        }
                        else{
                                w=exp(-lam1[n_i]*y1[n_i])
                                v=sample(uv[uv[,1]<=w,2],1,replace=T)
                                r_new=-log(v)/(lam2[n_i])
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
                }
                s3=1-s2
                s3[length(s2)]=(t1[n_i]<=c1[n_i])
                s11=(t1[n_i]<=c1[n_i])
                main[[n_i]]=cbind(id=n_i,y2=t2,s1=s11,s2,s3,x1=x1[n_i],x2=x2[n_i],y1=y1[n_i])
        }
        ##########################################

        main=data.frame(do.call("rbind",main))


        b1_sample=rep(b1,n_sample)
        b2_sample=rep(b2,n_sample)
        b3_sample=rep(b3,n_sample)
        b4_sample=rep(b4,n_sample)
        b01_sample=rep(b01,n_sample)
        b02_sample=rep(b02,n_sample)
        sigma_sample=rep(1,n_sample)
        theta_sample=rep(mu,n_sample)


        scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1) ###step size of MCMC###


        phiN_est_list=MH_joint_model(b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
                                     b01=as.double(b01_sample), b02=as.double(b02_sample),
                                     phi=as.double(phi),sigma=as.double(sigma_sample),
                                     theta=as.double(theta_sample), x1=as.double(main$x1), x2=as.double(main$x2),
                                     y1=as.double(main$y1),y2=as.double(main$y2), s1=as.integer(main$s1),
                                     s2=as.integer(main$s2), id=as.integer(main$id-1), no_pat_p=as.integer(length(unique(main$id))),
                                     n_sample=as.integer(n_sample),  N_p=as.integer(length(main$y1)),scale=as.double(scale))

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
        theta_sample=rep(theta_est1,n_sample)
        sigma_epsilon_sample=rep(sigma_epsilon,n_sample)
        mu_sample=rep(theta_est1,n_sample)

        scale=c(0.1,0.1,0.1,0.1,0.1,0.1,0.01,0.1,0.01)


        result=MH_dep(b1=as.double(b1_sample), b2=as.double(b2_sample),b3=as.double(b3_sample),b4=as.double(b4_sample),
                      b01=as.double(b01_sample), b02=as.double(b02_sample),
                      phi=as.double(phi),sigma=as.double(sigma_sample),mu=as.double(mu_sample),
                      theta=as.double(theta_sample),
                      rho=as.double(rho_sample), sigma_epsilon=as.double(sigma_epsilon_sample), x1=as.double(main$x1), x2=as.double(main$x2),
                      y1=as.double(main$y1),y2=as.double(main$y2), s1=as.integer(main$s1),
                      s2=as.integer(main$s2), id=as.integer(main$id-1), no_pat_p=as.integer(length(unique(main$id))),
                      n_sample=as.integer(n_sample),  N_p=as.integer(length(main$y1)),scale=as.double(scale),ni=as.integer(main$ni),maxni=as.integer(main$maxni))


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


        return(list(b01_est1,b02_est1,b1_est1,b2_est1,b3_est1,b4_est1,
                    b01_est,b02_est,b1_est,b2_est,b3_est,b4_est))
}

