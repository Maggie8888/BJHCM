#include <string>
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
//' Return paramter estimation in TJCM
//'

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
//' @export

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


double LL (double b1, double b2,double b3,double b4, double b01, double b02, NumericVector phi,
           double theta, NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,
           IntegerVector s1,IntegerVector s2,IntegerVector id, int N)
{
        double sum=0.0,lam1,lam2,U,V,u,v;
        for (int i = 0; i < N; i++){
                lam1=exp(b01+phi[id[i]]+b1*x1[i]+b2*x2[i]);
                lam2=exp(b02+phi[id[i]]+b3*x1[i]+b4*x2[i]);
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

double LL_id(double b1, double b2,double b3,double b4, double b01, double b02,
             double phi, double theta,NumericVector x1,NumericVector x2,NumericVector y1,
             NumericVector y2,IntegerVector s1,IntegerVector s2, int N)
{
        double sum=0.0,lam1,lam2,U,V,u,v;
        for (int i = 0; i < N; i++){
                lam1=exp(b01+phi+b1*x1[i]+b2*x2[i]);
                lam2=exp(b02+phi+b3*x1[i]+b4*x2[i]);
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

double target_phi(double b1, double b2, double b3, double b4, double b01, double b02,
                  NumericVector phi,double sigma,double theta,NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,
                  IntegerVector s1,IntegerVector s2,IntegerVector id,int no_patient,int N, double x1_new, int k)
{
        double sum = 0.0, func=0.0;
        int m;
        NumericVector temp=phi;
        for(m=0; m<no_patient; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=phi[m];
                }
        }
        func=func+sigma*temp[k]*temp[k];
        sum=LL(b1,b2,b3,b4,b01,b02,temp,theta,x1,x2, y1,y2, s1, s2,id,N)-0.5*func;
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

double metroplis_phi(double b1,double b2,double b3, double b4, double b01, double b02,NumericVector phi_curr,double sigma, double theta,
                     NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,
                     IntegerVector id,int no_patient, int N, double x1_prev, int j, double scale)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_phi(b1,b2,b3,b4,b01,b02,phi_curr,sigma,theta,x1,x2,y1,y2,s1,s2,id,no_patient,N,x1_curr,j);
        temp2=target_phi(b1,b2,b3,b4,b01,b02,phi_curr,sigma,theta,x1,x2,y1,y2,s1,s2,id,no_patient,N,x1_prev,j);
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

double metroplis_phi_id(double b1,double b2,double b3, double b4, double b01, double b02, double theta, double sigma,
                        NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,IntegerVector s1, IntegerVector s2,int N, double x1_prev, double scale)
{
        double lgratio=0, u ,temp1=0, temp2=0;
        double x1_curr, phi_curr=0.0;
        u=R::runif(0,1);
        x1_curr=x1_prev+scale*R::rnorm(0,1);
        temp1=LL_id(b1,b2,b3,b4,b01,b02,x1_curr,theta,x1,x2,y1,y2,s1,s2,N)+R::dnorm(x1_curr,0,sigma,1);
        temp2=LL_id(b1,b2,b3,b4,b01,b02,x1_prev,theta,x1,x2,y1,y2,s1,s2,N)+R::dnorm(x1_prev,0,sigma,1);
        lgratio=temp1-temp2;

        if (lgratio>=log(u))
        {
                phi_curr=x1_curr;
        }
        else
        {
                phi_curr=x1_prev;
        }
        return(phi_curr);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH_joint_model(NumericVector b1,NumericVector b2,NumericVector b3,NumericVector b4,NumericVector b01,
                    NumericVector b02,NumericVector phi,NumericVector sigma,NumericVector theta,
                    NumericVector x1,NumericVector x2,
                    NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,IntegerVector id,int no_pat_p,
                    int n_sample,int N_p,NumericVector scale){
        int i,j;
        double u;
        double y_new,lgratio,temp1,temp2;
        double sig_pat,x1_prev;
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
                y_new=b1[i-1]+scale[0]*R::rnorm(0,1);
                temp1=LL(y_new,b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i-1],b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b1[i]=y_new ;
                else b1[i]=b1[i-1];

                u=R::runif(0,1);
                y_new=b2[i-1]+scale[1]*R::rnorm(0,1);
                temp1=LL(b1[i],y_new,b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i],b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b2[i]=y_new ;
                else b2[i]=b2[i-1];

                u=R::runif(0,1);
                y_new=b3[i-1]+scale[2]*R::rnorm(0,1);
                temp1=LL(b1[i],b2[i],y_new,b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i],b2[i],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b3[i]=y_new ;
                else b3[i]=b3[i-1];

                u=R::runif(0,1);
                y_new=b4[i-1]+scale[3]*R::rnorm(0,1);
                temp1=LL(b1[i],b2[i],b3[i],y_new,b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i],b2[i],b3[i],b4[i-1],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b4[i]=y_new ;
                else b4[i]=b4[i-1];

                u=R::runif(0,1);
                y_new=b01[i-1]+scale[4]*R::rnorm(0,1);
                temp1=LL(b1[i],b2[i],b3[i],b4[i],y_new,b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i],b2[i],b3[i],b4[i],b01[i-1],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b01[i]=y_new ;
                else b01[i]=b01[i-1];


                u=R::runif(0,1);
                y_new=b02[i-1]+scale[5]*R::rnorm(0,1);
                temp1=LL(b1[i],b2[i],b3[i],b4[i],b01[i],y_new,phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i-1],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b02[i]=y_new ;
                else b02[i]=b02[i-1];


                /*********theta*******/
                u=R::runif(0,1);
                y_new=theta[i-1]+scale[6]*R::rnorm(0,1);
                if(y_new>0){
                        temp1=LL(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,y_new,x1, x2, y1,y2, s1, s2,id,N_p)+R::dgamma(y_new,0.001,1000,1);
                        temp2=LL(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,theta[i-1],x1, x2, y1,y2, s1, s2,id,N_p)+R::dgamma(theta[i-1],0.001,1000,1);
                        lgratio=temp1-temp2;

                        if (lgratio>=log(u))  theta[i]=y_new ;
                        else theta[i]=theta[i-1];
                }
                else theta[i]=theta[i-1];


                /************************/


                for (j = 0; j < no_pat_p; j++)
                {
                        x1_prev=phi[j];
                        metroplis_phi(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,sigma[i-1],theta[i],x1,x2,y1,y2,s1,s2,id,no_pat_p,N_p,x1_prev,j,scale[7]);
                        phi=center(phi,no_pat_p);

                }




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

                sig_pat=target_sigma(phi,no_pat_p);
                sigma[i]=R::rgamma(0.001+no_pat_p/2,1/(0.001+0.5*sig_pat));

        }


        PutRNGstate();
        return List::create(
                _["phi"]=phi,
                _["theta"]= theta,
                _["b1"]= b1,
                _["b2"]= b2,
                _["b3"]= b3,
                _["b4"]= b4,
                _["b01"]= b01,
                _["b02"]= b02,
                _["sigma"]= sigma);

}


