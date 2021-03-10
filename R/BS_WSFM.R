#' ####################################################################
#' ########### This function  return figure of BS #############
#' ############ input #########
#' idx: index control simulation settings

#' ##########  Output  ############
#' figure for BS
#'####################################################################
#'##########    usage ##############################

#' result<-readRDS("200_WSF_72.rds")
#' BS(72)

######################################################
######   generate data and fit the model   ##########
######################################################
#'
#' @param idx
#' @return graph for BS
#' @export
BS_figure<-function(idx){
  ##idx=1 for the scenario with sample size n=100 and theta=1,rho=0.1,death rate=0.6; idx=2 for the scenario with sample size n=100 and theta=1,death rate=0.5,rho=0.1 ##
  ##idx=3 for the scenario with sample size n=100 and theta=1,death rate=0.4 ,rho=0.1; idx=4 for the scenario with sample size n=100 and theta=1.5,death rate=0.6,rho=0.1 ##
  ##idx=5 for the scenario with sample size n=100 and theta=1.5 ,death rate=0.5,rho=0.1 ; idx=6 for the scenario with sample size n=100 and theta=1.5,death rate=0.4,rho=0.1 ##
  ##idx=7 for the scenario with sample size n=100 and theta=2 ,death rate=0.6,rho=0.1 ; idx=8 for the scenario with sample size n=100 and theta=2,death rate=0.5,rho=0.1 ##
  ##idx=9 for the scenario with sample size n=100 and theta=2,death rate=0.4,rho=0.1 ;

  ## idx=10 for the scenario with sample size n=200 and theta=1,death rate=0.6,rho=0.1
  ##idx=11 for the scenario with sample size n=200 and theta=1 ,death rate=0.5,rho=0.1 ; idx=12 for the scenario with sample size n=200 and theta=1,death rate=0.4,rho=0.1 ##
  ##idx=13 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.6,rho=0.1 ; idx=14 for the scenario with sample size n=200 and theta=1.5,death rate=0.5,rho=0.1 ##
  ##idx=15 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.4,rho=0.1 ; idx=16 for the scenario with sample size n=200 and theta=2,death rate=0.6,rho=0.1 ##
  ##idx=17 for the scenario with sample size n=200 and theta=2 ,death rate=0.5,rho=0.1 ; idx=18 for the scenario with sample size n=200 and theta=2,death rate=0.4,rho=0.1 ##

  ## idx=19 for the scenario with sample size n=400 and theta=1,death rate=0.6,rho=0.1
  ##idx=20 for the scenario with sample size n=400 and theta=1 ,death rate=0.5,rho=0.1 ; idx=21 for the scenario with sample size n=400 and theta=1,death rate=0.4,rho=0.1 ##
  ##idx=22 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.6,rho=0.1 ; idx=23 for the scenario with sample size n=400 and theta=1.5,death rate=0.5,rho=0.1 ##
  ##idx=24 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.4,rho=0.1 ; idx=25 for the scenario with sample size n=400 and theta=2,death rate=0.6,rho=0.1 ##
  ##idx=26 for the scenario with sample size n=400 and theta=2 ,death rate=0.5,rho=0.1 ; idx=27 for the scenario with sample size n=400 and theta=2,death rate=0.4,rho=0.1 ##


  ##idx=28 for the scenario with sample size n=100 and theta=1,death rate=0.6,rho=0.3; idx=29 for the scenario with sample size n=100 and theta=1,death rate=0.5,rho=0.3 ##
  ##idx=30 for the scenario with sample size n=100 and theta=1,death rate=0.4,rho=0.3; idx=31 for the scenario with sample size n=100 and theta=1.5,death rate=0.6,rho=0.3 ##
  ##idx=32 for the scenario with sample size n=100 and theta=1.5 ,death rate=0.5,rho=0.3 ; idx=33 for the scenario with sample size n=100 and theta=1.5,death rate=0.4,rho=0.3 ##
  ##idx=34 for the scenario with sample size n=100 and theta=2 ,death rate=0.6,rho=0.3 ; idx=35 for the scenario with sample size n=100 and theta=2,death rate=0.5,rho=0.3 ##
  ##idx=36 for the scenario with sample size n=100 and theta=2,death rate=0.4,rho=0.3 ;

  ##idx=37 for the scenario with sample size n=200 and theta=1,death rate=0.6,rho=0.3 ##
  ##idx=38 for the scenario with sample size n=200 and theta=1 ,death rate=0.5,rho=0.3 ;
  ##idx=39 for the scenario with sample size n=200 and theta=1,death rate=0.4,rho=0.3 ##
  ##idx=40 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.6,rho=0.3 ; idx=41 for the scenario with sample size n=200 and theta=1.5,death rate=0.5,rho=0.3 ##
  ##idx=42 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.4,rho=0.3 ; idx=43 for the scenario with sample size n=200 and theta=2,death rate=0.6,rho=0.3 ##
  ##idx=44 for the scenario with sample size n=200 and theta=2 ,death rate=0.5,rho=0.3 ; idx=45 for the scenario with sample size n=200 and theta=2,death rate=0.4,rho=0.3 ##

  ##idx=46 for the scenario with sample size n=400 and theta=1,death rate=0.6,rho=0.3 ##
  ##idx=47 for the scenario with sample size n=400 and theta=1 ,death rate=0.5,rho=0.3 ;
  ##idx=48 for the scenario with sample size n=400 and theta=1,death rate=0.4,rho=0.3 ##
  ##idx=49 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.6,rho=0.3 ; idx=50 for the scenario with sample size n=200 and theta=1.5,death rate=0.5,rho=0.3 ##
  ##idx=51 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.4,rho=0.3 ; idx=52 for the scenario with sample size n=200 and theta=2,death rate=0.6,rho=0.3 ##
  ##idx=53 for the scenario with sample size n=400 and theta=2 ,death rate=0.5,rho=0.3 ; idx=54  for the scenario with sample size n=200 and theta=2,death rate=0.4,rho=0.3 ##


  ##idx=55 for the scenario with sample size n=100 and theta=1,death rate=0.6,rho=0.5; idx=56 for the scenario with sample size n=100 and theta=1,death rate=0.5,rho=0.5 ##
  ##idx=57 for the scenario with sample size n=100 and theta=1,death rate=0.4 ,rho=0.5; idx=58 for the scenario with sample size n=100 and theta=1.5,death rate=0.6,rho=0.5 ##
  ##idx=59 for the scenario with sample size n=100 and theta=1.5 ,death rate=0.5,rho=0.5 ; idx=60 for the scenario with sample size n=100 and theta=1.5,death rate=0.4,rho=0.5 ##
  ##idx=61 for the scenario with sample size n=100 and theta=2 ,death rate=0.6,rho=0.5 ; idx=62 for the scenario with sample size n=100 and theta=2,death rate=0.5,rho=0.5 ##
  ##idx=63 for the scenario with sample size n=100 and theta=2,death rate=0.4,rho=0.5 ;

  ##idx=64 for the scenario with sample size n=200 and theta=1,death rate=0.6,rho=0.5 ##
  ##idx=65 for the scenario with sample size n=200 and theta=1 ,death rate=0.5,rho=0.5 ;
  ##idx=66 for the scenario with sample size n=200 and theta=1,death rate=0.4,rho=0.5 ##
  ##idx=67 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.6,rho=0.5 ; idx=68 for the scenario with sample size n=200 and theta=1.5,death rate=0.5,rho=0.5 ##
  ##idx=69 for the scenario with sample size n=200 and theta=1.5 ,death rate=0.4,rho=0.5 ; idx=70 for the scenario with sample size n=200 and theta=2,death rate=0.6,rho=0.5 ##
  ##idx=71 for the scenario with sample size n=200 and theta=2 ,death rate=0.5,rho=0.5 ; idx=72 for the scenario with sample size n=200 and theta=2,death rate=0.4,rho=0.5 ##

  ##idx=73 for the scenario with sample size n=400 and theta=1,death rate=0.6,rho=0.5 ##
  ##idx=74 for the scenario with sample size n=400 and theta=1 ,death rate=0.5,rho=0.5 ;
  ##idx=75 for the scenario with sample size n=400 and theta=1,death rate=0.4,rho=0.5 ##
  ##idx=76 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.6,rho=0.5 ; idx=77 for the scenario with sample size n=400 and theta=1.5,death rate=0.5,rho=0.5 ##
  ##idx=78 for the scenario with sample size n=400 and theta=1.5 ,death rate=0.4,rho=0.5 ; idx=79 for the scenario with sample size n=400 and theta=2,death rate=0.6,rho=0.5 ##
  ##idx=80 for the scenario with sample size n=400 and theta=2 ,death rate=0.5,rho=0.5 ; idx=81  for the scenario with sample size n=400 and theta=2,death rate=0.4,rho=0.5 ##


  bs_t=list()
  bs_t_terminal=list()
  bs_t_terminal_SF=list()

  bs_t_up=list()
  bs_t_terminal_up=list()
  bs_t_terminal_SF_up=list()
  bs_t_sd=list()
  bs_t_terminal_sd=list()
  bs_t_terminal_SF_sd=list()





  bs_t[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_sd[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t[[1]][,2]<-c(rep(0,8))
  bs_t[[1]][,1]<-c(rep(0,8))
  bs_t[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_sd[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t[[2]][,2]<-c(rep(0,5))
  bs_t[[2]][,1]<-c(rep(0,5))
  bs_t[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_sd[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t[[3]][,2]<-c(rep(0,2))
  bs_t[[3]][,1]<-c(rep(0,2))
  bs_t_terminal[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_terminal_sd[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_terminal[[1]][,2]<-c(rep(0,8))
  bs_t_terminal[[1]][,1]<-c(rep(0,8))
  bs_t_terminal[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_terminal_sd[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_terminal[[2]][,2]<-c(rep(0,5))
  bs_t_terminal[[2]][,1]<-c(rep(0,5))
  bs_t_terminal[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_terminal_sd[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_terminal[[3]][,2]<-c(rep(0,2))
  bs_t_terminal[[3]][,1]<-c(rep(0,2))
  bs_t_terminal_SF[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_terminal_SF_sd[[1]]<-matrix(c(rep(0,16)),nrow=8,ncol=2)
  bs_t_terminal_SF[[1]][,2]<-c(rep(0,8))
  bs_t_terminal_SF[[1]][,1]<-c(rep(0,8))
  bs_t_terminal_SF[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_terminal_SF_sd[[2]]<-matrix(c(rep(0,10)),nrow=5,ncol=2)
  bs_t_terminal_SF[[2]][,2]<-c(rep(0,5))
  bs_t_terminal_SF[[2]][,1]<-c(rep(0,5))
  bs_t_terminal_SF[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_terminal_SF_sd[[3]]<-matrix(c(rep(0,4)),nrow=2,ncol=2)
  bs_t_terminal_SF[[3]][,2]<-c(rep(0,2))
  bs_t_terminal_SF[[3]][,1]<-c(rep(0,2))




  bs_t_up[[1]]<-matrix(c(rep(0,8*length(result))),nrow=length(result),ncol=8)
  bs_t_up[[2]]<-matrix(c(rep(0,5*length(result))),nrow=length(result),ncol=5)
  bs_t_up[[3]]<-matrix(c(rep(0,2*length(result))),nrow=length(result),ncol=2)

  bs_t_terminal_up[[1]]<-matrix(c(rep(0,8*length(result))),nrow=length(result),ncol=8)
  bs_t_terminal_up[[2]]<-matrix(c(rep(0,5*length(result))),nrow=length(result),ncol=5)
  bs_t_terminal_up[[3]]<-matrix(c(rep(0,2*length(result))),nrow=length(result),ncol=2)

  bs_t_terminal_SF_up[[1]]<-matrix(c(rep(0,8*length(result))),nrow=length(result),ncol=8)
  bs_t_terminal_SF_up[[2]]<-matrix(c(rep(0,5*length(result))),nrow=length(result),ncol=5)
  bs_t_terminal_SF_up[[3]]<-matrix(c(rep(0,2*length(result))),nrow=length(result),ncol=2)


  for (i in 1:length(result)){
    for (j in 1:8) {
      temp1<-ifelse(is.null(result[[i]][[1]][[1]][j,2]),0,result[[i]][[1]][[1]][j,2])
      bs_t_up[[1]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[1]][j,2]),0,result[[i]][[2]][[1]][j,2])
      bs_t_terminal_up[[1]][i,j]<-temp2
    }
  }


  for (i in 1:length(result)){
    for (j in 1:5) {
      temp1<-ifelse(is.null(result[[i]][[1]][[2]][j,2]),0,result[[i]][[1]][[2]][j,2])
      bs_t_up[[2]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[2]][j,2]),0,result[[i]][[2]][[2]][j,2])
      bs_t_terminal_up[[2]][i,j]<-temp2
    }
  }
  for (i in 1:length(result)){
    for (j in 1:2) {
      temp1<-ifelse(is.null(result[[i]][[1]][[3]][j,2]),0,result[[i]][[1]][[3]][j,2])
      bs_t_up[[3]][i,j]<-temp1
      temp2<-ifelse(is.null(result[[i]][[2]][[3]][j,2]),0,result[[i]][[2]][[3]][j,2])
      bs_t_terminal_up[[3]][i,j]<-temp2
    }
  }

  bs_t_sd[[1]][1:8,2]<-apply(bs_t_up[[1]], 2, sd)
  bs_t_terminal_sd[[1]][1:8,2]<-apply(bs_t_terminal_up[[1]], 2, sd)
  bs_t_sd[[2]][1:5,2]<-apply(bs_t_up[[2]], 2, sd)
  bs_t_terminal_sd[[2]][1:5,2]<-apply(bs_t_terminal_up[[2]], 2, sd)
  bs_t_sd[[3]][1:2,2]<-apply(bs_t_up[[3]], 2, sd)
  bs_t_terminal_sd[[3]][1:2,2]<-apply(bs_t_terminal_up[[3]], 2, sd)

  bs_t[[1]][1:8,2]<-apply(bs_t_up[[1]], 2, mean)
  bs_t_terminal[[1]][1:8,2]<-apply(bs_t_terminal_up[[1]], 2, mean)
  bs_t[[2]][1:5,2]<-apply(bs_t_up[[2]], 2, mean)
  bs_t_terminal[[2]][1:5,2]<-apply(bs_t_terminal_up[[2]], 2, mean)
  bs_t[[3]][1:2,2]<-apply(bs_t_up[[3]], 2, mean)
  bs_t_terminal[[3]][1:2,2]<-apply(bs_t_terminal_up[[3]], 2, mean)


  bs_t[[1]][1:8,1]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t[[2]][1:5,1]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t[[3]][1:2,1]<-c(0.1,0.09)
  bs_t_sd[[1]][1:8,1]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t_sd[[2]][1:5,1]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t_sd[[3]][1:2,1]<-c(0.1,0.09)


  bs_t_terminal[[1]][1:8]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t_terminal[[2]][1:5]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t_terminal[[3]][1:2]<-c(0.1,0.09)
  bs_t_terminal_sd[[1]][1:8]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t_terminal_sd[[2]][1:5]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t_terminal_sd[[3]][1:2]<-c(0.1,0.09)


  for (i in 1:length(result)){
    for (j in 1:8) {

      temp<-ifelse(is.null(result[[i]][[3]][[1]][j,2]),0,result[[i]][[3]][[1]][j,2])
      bs_t_terminal_SF_up[[1]][i,j]<-temp
    }
  }

  for (i in 1:length(result)){
    for (j in 1:5) {
      temp<-ifelse(is.null(result[[i]][[3]][[2]][j,2]),0,result[[i]][[3]][[2]][j,2])
      bs_t_terminal_SF_up[[2]][i,j]<-temp
    }
  }
  for (i in 1:length(result)){
    for (j in 1:2) {
      temp<-ifelse(is.null(result[[i]][[3]][[3]][j,2]),0,result[[i]][[3]][[3]][j,2])
      bs_t_terminal_SF_up[[3]][i,j]<-temp
    }
  }

  bs_t_sd[[1]][1:8,2]<-apply(bs_t_up[[1]], 2, sd)
  bs_t_terminal_sd[[1]][1:8,2]<-apply(bs_t_terminal_up[[1]], 2, sd)
  bs_t_terminal_SF_sd[[1]][8:1,2]<-apply(bs_t_terminal_SF_up[[1]], 2, sd)
  bs_t_sd[[2]][1:5,2]<-apply(bs_t_up[[2]], 2, sd)
  bs_t_terminal_sd[[2]][1:5,2]<-apply(bs_t_terminal_up[[2]], 2, sd)
  bs_t_terminal_SF_sd[[2]][5:1,2]<-apply(bs_t_terminal_SF_up[[2]], 2, sd)
  bs_t_sd[[3]][1:2,2]<-apply(bs_t_up[[3]], 2, sd)
  bs_t_terminal_sd[[3]][1:2,2]<-apply(bs_t_terminal_up[[3]], 2, sd)
  bs_t_terminal_SF_sd[[3]][2:1,2]<-apply(bs_t_terminal_SF_up[[3]], 2, sd)


  bs_t_terminal_SF[[1]][1:8,2]<-apply(bs_t_terminal_SF_up[[1]], 2, mean)
  bs_t_terminal_SF[[2]][1:5,2]<-apply(bs_t_terminal_SF_up[[2]], 2, mean)
  bs_t_terminal_SF[[3]][1:2,2]<-apply(bs_t_terminal_SF_up[[3]], 2, mean)
  bs_t_terminal_SF[[1]][1:8]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t_terminal_SF[[2]][1:5]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t_terminal_SF[[3]][1:2]<-c(0.1,0.09)
  bs_t_terminal_SF_sd[[1]][1:8]<-c(0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03)
  bs_t_terminal_SF_sd[[2]][1:5]<-c(0.1,0.09,0.08,0.07,0.06)
  bs_t_terminal_SF_sd[[3]][1:2]<-c(0.1,0.09)

  setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/cluster/")
  #setwd("/Users/mengluliang/Dropbox/biostatisitcs_PSU/RA/2020Spring/BayesianJointModel/03preresults2/")
  res<-list(bs_t,bs_t_terminal,bs_t_terminal_SF,bs_t_sd,bs_t_terminal_sd,bs_t_terminal_SF_sd)
  #sink(paste0("output",idx,".txt"))
  #print(res)
  #sink()
  cat(capture.output(print(res), file=paste0("output_WSFM",idx,".txt")))
#cat(capture.output(print(bs_t), file=paste0("BS_terminal",idx,".txt")))
  #cat(capture.output(print(bs_t_terminal), file=paste0("BS_terminal_old",idx,".txt")))
  #cat(capture.output(print(bs_t_sd), file=paste0("BS_terminal_SD",idx,".txt")))
  #cat(capture.output(print(bs_t_terminal_sd), file=paste0("BS_terminal_SD_old",idx,".txt")))

  temp<-as.data.frame(cbind(bs_t[[1]],bs_t_sd[[1]][,2],bs_t_terminal[[1]][,2],
              bs_t_terminal_sd[[1]][,2],bs_t_terminal_SF[[1]][,2],
              bs_t_terminal_SF_sd[[1]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_0.03_WSFM",idx,".csv"))

  temp<-as.data.frame(cbind(bs_t[[2]],bs_t_sd[[2]][,2],bs_t_terminal[[2]][,2],
                            bs_t_terminal_sd[[2]][,2],bs_t_terminal_SF[[2]][,2],
                            bs_t_terminal_SF_sd[[2]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_0.06_WSFM",idx,".csv"))

  temp<-as.data.frame(cbind(bs_t[[3]],bs_t_sd[[3]][,2],bs_t_terminal[[3]][,2],
                            bs_t_terminal_sd[[3]][,2],bs_t_terminal_SF[[3]][,2],
                            bs_t_terminal_SF_sd[[3]][,2]))
  colnames(temp)<-c("t","terminal_mean","terminal_sd","terminal_mean_old","terminal_sd_old")
  write.csv(temp, file=paste0("BS_terminal_0.09_WSFM",idx,".csv"))


  png(paste0("BS_WSFM_",idx,".png"),width=900, height=900)

  t_s=c(0.03,0.06,0.09)

  for(t_prime in t_s){

    idxx=which(t_s==t_prime)


    if(t_prime==t_s[1]) {

      plot(bs_t[[1]][1:8],bs_t[[1]][,2],xlab="",type="o",xlim=c(0,0.1),ylim=c(0,0.35),main="",lty=1,col = "red", ylab = "",pch=16,lwd=2)
      #arrows(x0=bs_t[[1]][1:8],y0=bs_t[[1]][,2]-bs_t_sd[[1]][,2], x1=bs_t[[1]][1:8], y1=bs_t[[1]][,2]+bs_t_sd[[1]][,2], code=3, angle=90, length=0.1,col = "red")
      legend("topleft",c("t'=0.03","t'=0.06","t'=0.09"),lty=1:3,lwd=1,pch=1,col = "red",title= "Time-varying JFCM",cex = 2)
      lines(bs_t_terminal[[1]],type="o",lty=1,lwd=2)
      lines(bs_t_terminal_SF[[1]],type="o",lty=1,col="blue",lwd=2)
      mtext(side=1, line=2, "time", font=3,cex=1.2)
      mtext(side=2, line=3, "Brier score",font=3,cex=1.2)
      #arrows(x0=bs_t_terminal[[1]][1:8],y0=bs_t_terminal[[1]][,2]-bs_t_terminal_sd[[1]][,2], x1=bs_t_terminal[[1]][1:8], y1=bs_t_terminal[[1]][,2]+bs_t_terminal_sd[[1]][,2], code=3, angle=90, length=0.1)
    }
    else {
      lines(bs_t[[2]],type="o",lty=2,col = "red",lwd=2)
      #arrows(x0=bs_t[[2]][1:5],y0=bs_t[[2]][,2]-bs_t_sd[[2]][,2], x1=bs_t[[2]][1:5], y1=bs_t[[2]][,2]+bs_t_sd[[2]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t[[3]],type="o",lty=3,col = "red",lwd=2)
      #arrows(x0=bs_t[[3]][1:2],y0=bs_t[[3]][,2]-bs_t_sd[[3]][,2], x1=bs_t[[3]][1:2], y1=bs_t[[3]][,2]+bs_t_sd[[3]][,2], code=3, angle=90, length=0.1,col = "red")
      lines(bs_t_terminal[[2]],type="o",lty=2,lwd=2)
      #arrows(x0=bs_t_terminal[[2]][1:5],y0=bs_t_terminal[[2]][,2]-bs_t_terminal_sd[[2]][,2], x1=bs_t_terminal[[2]][1:5], y1=bs_t_terminal[[2]][,2]+bs_t_terminal_sd[[2]][,2], code=3, angle=90, length=0.1)
      lines(bs_t_terminal[[3]],type="o",lty=3,lwd=2)
      #arrows(x0=bs_t_terminal[[3]][1:2],y0=bs_t_terminal[[3]][,2]-bs_t_terminal_sd[[3]][,2], x1=bs_t_terminal[[3]][1:2], y1=bs_t_terminal[[3]][,2]+bs_t_terminal_sd[[3]][,2], code=3, angle=90, length=0.1)
      lines(bs_t_terminal_SF[[2]],type="o",lty=2,col="blue",lwd=2)
      lines(bs_t_terminal_SF[[3]],type="o",lty=3,col="blue",lwd=2)
      }
    abline(v=t_prime,lty=2)
    if(t_prime==t_s[1]) {axis(side=1,seq(0.01,0.10,0.01))
      axis(side=2, seq(0,0.35,0.05))}
  }

  legend("bottomleft",c("t'=0.03","t'=0.06","t'=0.09"),lty=1:3,lwd=1,pch=1,title= "JFCM",cex = 2)
  legend("left",c("t'=0.03","t'=0.06","t'=0.09"),lty=1:3,lwd=1,pch=1,title= "SFM",cex = 2,col="blue")
  ###################################################################################

  dev.off()

  return(list(bs_t, bs_t_terminal,bs_t_terminal_SF,bs_t_sd,bs_t_terminal_sd,bs_t_terminal_SF_sd))
  ###################################################################################

}


