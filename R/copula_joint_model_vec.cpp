#include <string>
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <list> 
#include <iterator> 
//#include <bits/stdc++.h> 
#include <vector> 
#include <random>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <list> 
#include <stdlib.h>
using namespace std; 
using namespace arma;
using namespace Rcpp;
#include <iostream>

//using std::swap;
//using std::cout;
//using std::endl;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

NumericVector center(NumericVector x,int n_region){
  double sum=0.0, avg;
  for (int i=0;i<n_region;i++){
    sum=sum+x[i];
  }
  avg=sum/n_region;
  for (int j=0;j<n_region;j++){
    x[j]=x[j]-avg;
  }
  return x;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double LL (NumericVector b1,NumericVector b2, double b01, double b02, NumericVector phi, 
           double theta, NumericMatrix x,NumericVector y1,NumericVector y2,
           IntegerVector s1,IntegerVector s2,IntegerVector id, int N, int p,
           IntegerVector ni,IntegerVector maxni)
{
  double sum=0.0,lam1=0.0,lam2=0.0,U,V,u,v;
  NumericVector lam1_vec(p),lam2_vec(p);
  for (int i = 0; i < N; i++){
    for (int j=0; j<p; j++){
      
      lam1_vec[j]=b01+phi[id[i]]+b1[j]*x(i,j);
      lam2_vec[j]=b02+phi[id[i]]+b2[j]*x(i,j);
      lam1=lam1+lam1_vec[j];
      lam2=lam2+lam2_vec[j];
    }
    lam1=exp(lam1);
    lam2=exp(lam2);
    U=exp(-lam1*y1[i]);
    V=exp(-lam2*y2[i]);
    u=lam1*exp(-lam1*y1[i]);
    v=lam2*exp(-lam2*y2[i]);
    
    if(s1[i]==1&&s2[i]==1&& ni[i]!=maxni[i]){
      sum=sum+log(pow(U,-theta-1)*pow(V,-theta-1)*(1+theta))+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v);
    }
    else if (s1[i]==1&&s2[i]==1&& ni[i]==maxni[i]){
      sum=sum+log(pow(U,-theta-1)*pow(V,-theta-1)*(1+theta))+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v)+log(u);
    }
    else if(s1[i]==0&&s2[i]==0){
      sum=sum-1/theta*log(pow(U,-theta)+pow(V,-theta)-1);
    }
    else if(s1[i]==1&&s2[i]==0){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(U)+log(u);
    }
    else if(s1[i]==0&&s2[i]==1&&ni[i]!=maxni[i]){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(V)+log(v)-log(U);
    }
    else if(s1[i]==0&&s2[i]==1&& ni[i]==maxni[i]){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(V)+log(v);
    }
  }
  return sum;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


double LLi(NumericVector b1,NumericVector b2, double b01, double b02, NumericVector phi, 
           double theta, NumericMatrix x,NumericVector y1,NumericVector y2,
           IntegerVector s1,IntegerVector s2,IntegerVector id, int N, int p,
           IntegerVector j_start,IntegerVector j_end,IntegerVector ni,IntegerVector maxni)
{
  double sum=0.0,lam1=0.0,lam2=0.0,U,V,u,v;
  NumericVector lam1_vec(p),lam2_vec(p);
  int j_start_tmp=j_start[id[0]];
  int j_end_tmp=j_start[id[0]];
  for (int i = (j_start_tmp-1); i < j_end_tmp; i++){
    for (int j=0; j<p; j++){
      
      lam1_vec[j]=b01+phi[id[i]]+b1[j]*x(i,j);
      lam2_vec[j]=b02+phi[id[i]]+b2[j]*x(i,j);
      lam1=lam1+lam1_vec[j];
      lam2=lam2+lam2_vec[j];
    }
    lam1=exp(lam1);
    lam2=exp(lam2);
    U=exp(-lam1*y1[i]);
    V=exp(-lam2*y2[i]);
    u=lam1*exp(-lam1*y1[i]);
    v=lam2*exp(-lam2*y2[i]);
    
    if(s1[i]==1&&s2[i]==1&& ni[i]!=maxni[i]){
      sum=sum+log(pow(U,-theta-1)*pow(V,-theta-1)*(1+theta))+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v);
    }
    else if (s1[i]==1&&s2[i]==1&& ni[i]==maxni[i]){
      sum=sum+log(pow(U,-theta-1)*pow(V,-theta-1)*(1+theta))+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v)+log(u);
    }
    else if(s1[i]==0&&s2[i]==0){
      sum=sum-1/theta*log(pow(U,-theta)+pow(V,-theta)-1);
    }
    else if(s1[i]==1&&s2[i]==0){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(U)+log(u);
    }
    else if(s1[i]==0&&s2[i]==1&&ni[i]!=maxni[i]){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(V)+log(v)-log(U);
    }
    else if(s1[i]==0&&s2[i]==1&& ni[i]==maxni[i]){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(V)+log(v);
    }
  }
  return sum;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LL_id(NumericVector b1,NumericVector b2, double b01, double b02, 
             double phi,double theta,NumericMatrix x,NumericVector y1,
             NumericVector y2,IntegerVector s1,IntegerVector s2, int N, int p)
{
  double sum=0.0,lam1=0.0,lam2=0.0,U,V,u,v;
  NumericVector lam1_vec(p),lam2_vec(p);
  for (int i = 0; i < N; i++){
    for (int j=0; j<p; j++){
      
      lam1_vec[j]=b01+phi+b1[j]*x(i,j);
      lam2_vec[j]=b02+phi+b2[j]*x(i,j);
      lam1=lam1+lam1_vec[j];
      lam2=lam2+lam2_vec[j];
    }
    lam1=exp(lam1);
    lam2=exp(lam2);
    U=exp(-lam1*y1[i]);
    V=exp(-lam2*y2[i]);
    u=lam1*exp(-lam1*y1[i]);
    v=lam2*exp(-lam2*y2[i]);
    
    
    if(s1[i]==1&&s2[i]==1){
      sum=sum+log(pow(U,-theta-1)*pow(V,-theta-1)*(1+theta))+(-1/theta-2)*log(pow(U,-theta)+pow(V,-theta)-1)+log(v);
    }
    else if(s1[i]==0&&s2[i]==0){
      sum=sum-1/theta*log(pow(U,-theta)+pow(V,-theta)-1);
    }
    else if(s1[i]==1&&s2[i]==0){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(U)+log(u);
    }
    else if(s1[i]==0&&s2[i]==1){
      sum=sum+(-1/theta-1)*log(pow(U,-theta)+pow(V,-theta)-1)-(theta+1)*log(V)+log(v)-log(U);
    }
  }
  return sum;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_phi(NumericVector b1, NumericVector b2, double b01, double b02,
                  NumericVector phi,double sigma,double theta,NumericMatrix x,NumericVector y1,NumericVector y2,
                  IntegerVector s1,IntegerVector s2,IntegerVector id,int no_patient,int N, double x1_new, int k,int p,
                  IntegerVector j_start,IntegerVector j_end,IntegerVector ni,IntegerVector maxni)
{
  double sum = 0.0, func=0.0;
  int m;
  NumericVector tem(no_patient);
  tem=phi;
  for(m=0; m<no_patient; m++)
  {
    if(m==k) {
      tem[m]=x1_new;
    }
    else
    {
      tem[m]=phi[m];
    }
  }
  func=func+sigma*tem[k]*tem[k];
  sum=LLi(b1,b2,b01,b02,tem,theta,x,y1,y2, s1, s2,id,N,p,j_start,j_end,ni,maxni)-0.5*func;
  return (sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_sigma(NumericVector phi, int no_patient)
{
  double sum=0.0;
  for(int i=0; i<no_patient; i++)
  {
    sum=sum+phi[i]*phi[i];
  }
  return(sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double rand_gen() {
  // return a uniformly distributed random value
  return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double normalRandom() {
  // return a normally distributed random value
  double v1=rand_gen();
  double v2=rand_gen();
  return cos(2*3.14*v2)*sqrt(-2.*log(v1));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_phi(NumericVector b1,NumericVector b2,double b01, double b02,NumericVector phi_curr,double sigma, double theta,
                     NumericMatrix x,NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,
                     IntegerVector id,int no_patient, int N, double x1_prev, int j, double scale,int p,
                     IntegerVector j_start,IntegerVector j_end,IntegerVector ni,IntegerVector maxni)
{
  double lgratio=0,temp1=0,temp2=0;
  double x1_curr;
  double u=R::runif(0,1);
  x1_curr=x1_prev+scale*R::rnorm(0,1);
  temp1=target_phi(b1,b2,b01,b02,phi_curr,sigma,theta,x,y1,y2,s1,s2,id,no_patient,N,x1_curr,j,p,j_start,j_end,ni,maxni);
  temp2=target_phi(b1,b2,b01,b02,phi_curr,sigma,theta,x,y1,y2,s1,s2,id,no_patient,N,x1_prev,j,p,j_start,j_end,ni,maxni);
  lgratio=temp1-temp2;
  
  if (lgratio>=log(u))
  {
    phi_curr[j]=x1_curr;
  }
  else
  {
    phi_curr[j]=x1_prev;
  }
  
  return  phi_curr[j];
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_phi_id(NumericVector b1,NumericVector b2,double b01, double b02,double sigma, double theta,
                     NumericMatrix x,NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,
                     int N, double x1_prev,double scale,int p)
{
  double lgratio=0,temp1=0,temp2=0;
  double x1_curr,phi_curr=0.0;
  double u=R::runif(0,1);
  x1_curr=x1_prev+scale*R::rnorm(0,1);
  temp1=LL_id(b1,b2,b01,b02,x1_curr,theta,x,y1,y2,s1,s2,N,p)+R::dnorm(x1_curr,0,sigma,1);
  temp2=LL_id(b1,b2,b01,b02,x1_prev,theta,x,y1,y2,s1,s2,N,p)+R::dnorm(x1_prev,0,sigma,1);
  lgratio=temp1-temp2;
  
  if (lgratio>=log(u))
  {
    phi_curr=x1_curr;
  }
  else
  {
    phi_curr=x1_prev;
  }
  
  return  phi_curr;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH_joint_model(NumericMatrix b1,NumericMatrix b2,NumericVector b01, 
            NumericVector b02,NumericVector phi,NumericVector sigma,NumericVector theta,
            NumericMatrix x,NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,IntegerVector id,int no_pat_p, 
            int n_sample,int N_p,NumericVector scale, int p,
            IntegerVector j_start,IntegerVector j_end,IntegerVector ni,IntegerVector maxni){
  int i,j;
  double u,sig_pat,x1_prev,temp1,temp2,lgratio;
  NumericVector y_new1(p);
  double y_new2;
  GetRNGstate();
  printf("%d\n",n_sample);
  for (i = 1; i < n_sample; i++)
  {
    
    //double sd = 1.0;
    //double mean = 0.0;
    //double rnorm = normalRandom()*sd+mean;
    //y_new=b1[i-1]+scale[0]*rnorm;
    printf("%d\n",i);
    u=R::runif(0,1);

    for (int l=0;l<p;l++){
      y_new1[l]=b1(i-1,l)+scale[0]*R::rnorm(0,1);
    }
    
    temp1=LL(y_new1,b2.row(i-1),b01[i-1],b02[i-1],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    temp2=LL(b1.row(i-1),b2.row(i-1),b01[i-1],b02[i-1],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    lgratio=temp1-temp2;
    if (lgratio>=log(u))  b1.row(i)=y_new1;
    else b1.row(i)=b1.row(i-1);
    
    
    for (int l=0;l<p;l++){
      y_new1[l]=b2(i-1,l)+scale[1]*R::rnorm(0,1);
    }
    temp1=LL(b1.row(i),y_new1,b01[i-1],b02[i-1],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    temp2=LL(b1.row(i),b2.row(i-1),b01[i-1],b02[i-1],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    lgratio=temp1-temp2;
    if (lgratio>=log(u))  b2.row(i)=y_new1;
    else b2.row(i)=b2.row(i-1);
    
    
    y_new2=b01[i-1]+scale[4]*R::rnorm(0,1);
    temp1=LL(b1.row(i),b2.row(i),y_new2,b02[i-1],phi,theta[i-1],x, y1,y2,s1,s2,id,N_p,p,ni,maxni);
    temp2=LL(b1.row(i),b2.row(i),b01[i-1],b02[i-1],phi,theta[i-1],x,y1,y2,s1,s2,id,N_p,p,ni,maxni);
    lgratio=temp1-temp2;
    if (lgratio>=log(u))  b01[i]=y_new2 ;
    else b01[i]=b01[i-1];
    
  
    y_new2=b02[i-1]+scale[5]*R::rnorm(0,1);
    temp1=LL(b1.row(i),b2.row(i),b01[i],y_new2,phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    temp2=LL(b1.row(i),b2.row(i),b01[i],b02[i-1],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni);
    lgratio=temp1-temp2;
    if (lgratio>=log(u))  b02[i]=y_new2 ;
    else b02[i]=b02[i-1];
    
    
    /*********theta*******/
   
    y_new2=theta[i-1]+scale[6]*R::rnorm(0,1);
    if(y_new2>0){
      temp1=LL(b1.row(i),b2.row(i),b01[i],b02[i],phi,y_new2,x, y1,y2, s1, s2,id,N_p,p,ni,maxni)+R::dgamma(y_new2,0.001,1000,1);
      temp2=LL(b1.row(i),b2.row(i),b01[i],b02[i],phi,theta[i-1],x, y1,y2, s1, s2,id,N_p,p,ni,maxni)+R::dgamma(theta[i-1],0.001,1000,1);
      lgratio=temp1-temp2;
      
      if (lgratio>=log(u))  theta[i]=y_new2 ;
      else theta[i]=theta[i-1];
    }
    else theta[i]=theta[i-1];
    
    /************************/
    
    
    for (j = 0; j < no_pat_p; j++)
    {
     
      x1_prev=phi[j];
      metroplis_phi(b1.row(i),b2.row(i),b01[i],b02[i],phi,sigma[i-1],theta[i],x,y1,y2,s1,s2,id,no_pat_p,N_p,x1_prev,j,scale[7],p,
                    j_start,j_end,ni,maxni);
      
     phi=center(phi,no_pat_p); 
    }
     sig_pat=target_sigma(phi,no_pat_p);
    
    /*u=runif(0,1);
     y_new=rho[i-1]+scale[8]*rnorm(0,1);
     
     if(y_new>0&&y_new<1){
     temp1=log(sqrt(1-y_new*y_new))-0.5*target_sigma_epsilon(theta,N,ni,y_new)/sigma_epsilon[i-1];
     temp2=log(sqrt(1-rho[i-1]*rho[i-1]))-0.5*target_sigma_epsilon(theta,N,ni,rho[i-1])/sigma_epsilon[i-1];
     lgratio=temp1-temp2;
     if (lgratio>=log(u))  rho[i]=y_new ;
     else rho[i]=rho[i-1];
     
     }
     else rho[i]=rho[i-1];*/
    
    
   //sigma[i]=R::rgamma(0.001+no_pat_p/2,1/(0.001+0.5*sig_pat));
   sigma[i]=R::rgamma(0.001+no_pat_p/2,1/(0.001+0.5*sig_pat));
   //sigma[i]=R::rgamma(0.001+no_pat_p/2,1/(0.001+0.5*0.01));
  }
  
  
  PutRNGstate();
  return List::create(
    _["phi"]=phi,
    _["theta"]= theta,
    _["b1"]= b1,
    _["b2"]= b2,
    _["b01"]= b01,
    _["b02"]= b02,
    _["sigma"]= sigma);
  
}


