#include <string>
#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
//' Return parameters estimation in TVJHCM
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

NumericVector center_dep(NumericVector x,int n_region){
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

double LL_dep(double b1, double b2,double b3,double b4, double b01, double b02,NumericVector phi,
              double mu,NumericVector theta,NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,
              IntegerVector s1, IntegerVector s2,IntegerVector id, int N){
        double sum=0.0,lam1,lam2,U,V,u,v;
        for (int i = 0; i < N; i++){
                lam1=exp(b01+phi[id[i]]+b1*x1[i]+b2*x2[i]);
                lam2=exp(b02+phi[id[i]]+b3*x1[i]+b4*x2[i]);
                U=exp(-lam1*y1[i]);
                V=exp(-lam2*y2[i]);
                u=lam1*exp(-lam1*y1[i]);
                v=lam2*exp(-lam2*y2[i]);

                if(s1[i]==1 & s2[i]==1){
                        sum=sum+log(pow(U,-mu*exp(theta[i])-1)*pow(V,-mu*exp(theta[i])-1)*(1+mu*exp(theta[i])))+(-1/(mu*exp(theta[i]))-2)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)+log(v);
                }
                else if(s1[i]==0 & s2[i]==0){
                        sum=sum-1/(mu*exp(theta[i]))*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1);
                }
                else if(s1[i]==1 & s2[i]==0){
                        sum=sum+(-1/(mu*exp(theta[i]))-1)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)-(mu*exp(theta[i])+1)*log(U)+log(u);
                }
                else if(s1[i]==0 & s2[i]==1){
                        sum=sum+(-1/(mu*exp(theta[i]))-1)*log(pow(U,-mu*exp(theta[i]))+pow(V,-mu*exp(theta[i]))-1)-(mu*exp(theta[i])+1)*log(V)+log(v)-log(U);
                }
        }
        return sum;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double LL_id_dep(double b1, double b2,double b3,double b4, double b01, double b02, double phi, double sigma,
                 double mu,NumericVector theta,NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,
                 IntegerVector s1, IntegerVector s2, int N)
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
        return sum;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_phi_dep(double b1, double b2, double b3, double b4, double b01, double b02,
                      NumericVector phi, double sigma, double mu,NumericVector theta,NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,
                      IntegerVector s1,IntegerVector s2,IntegerVector id,int no_patient, int N, double x1_new, int k)
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
        func=func+temp[k]*temp[k]/sigma;
        sum=LL_dep(b1,b2,b3,b4,b01,b02,temp,mu,theta,x1, x2, y1,y2, s1, s2,id,N)-0.5*func;
        return (sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double target_sigma_dep(NumericVector phi, int no_patient)
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

double target_sigma_epsilon_dep(NumericVector theta, int N,IntegerVector ni, double rho)
{
        double sum=0.0;
        for(int i=1; i<N; i++)
        {
                if(ni[i]==1) sum=sum+(1-rho*rho)*theta[i]*theta[i];
                else sum=sum+(theta[i]-rho*theta[i-1])*(theta[i]-rho*theta[i-1]);
        }

        return(sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double theta_pow_sum_dep(NumericVector theta,int N,IntegerVector ni)
{
        double sum=0.0;
        for(int i=0; i<N; i++)
        {
                if(ni[i]>=2) sum=sum+theta[i]*theta[i];
        }
        return(sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double theta_lag_sum_dep(NumericVector theta,int N,IntegerVector ni)
{
        double sum=0.0;
        for(int i=1; i<N; i++)
        {
                if(ni[i]>=2) sum=sum+theta[i]*theta[i-1];

        }
        return(sum);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double rand_gen_dep() {
        // return a uniformly distributed random value
        return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double normalRandom_dep() {
        // return a normally distributed random value
        double v1=rand_gen_dep();
        double v2=rand_gen_dep();
        return cos(2*3.14*v2)*sqrt(-2.*log(v1));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_phi_dep(double b1,double b2,double b3, double b4, double b01, double b02,NumericVector phi_curr,double sigma, double mu,
                         NumericVector theta,NumericVector x1,NumericVector x2,NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,
                         IntegerVector id,int no_patient, int N, double x1_prev, int j, double scale)
{
        double lgratio=0,temp1=0,temp2=0;
        double x1_curr;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom_dep()*sd+mean;
        x1_curr=x1_prev+scale*rnorm;
        temp1=target_phi_dep(b1,b2,b3,b4,b01,b02,phi_curr,sigma,mu,theta,x1,x2,y1,y2,s1,s2,id,no_patient,N,x1_curr,j);
        temp2=target_phi_dep(b1,b2,b3,b4,b01,b02,phi_curr,sigma,mu,theta,x1,x2,y1,y2,s1,s2,id,no_patient,N,x1_prev,j);
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

double metroplis_phi_id_dep(double b1,double b2,double b3, double b4, double b01, double b02,
                            double sigma, double mu,NumericVector theta,NumericVector x1,
                            NumericVector x2,NumericVector y1,NumericVector y2,IntegerVector s1,
                            IntegerVector s2, int N,double x1_prev, double scale)
{
        double lgratio=0, u ,temp1=0, temp2=0;
        double x1_curr, phi_curr=0.0;
        u=R::runif(0,1);
        x1_curr=x1_prev+scale*R::rnorm(0,1);
        temp1=LL_id_dep(b1,b2,b3,b4,b01,b02,x1_curr,sigma,mu,theta,x1,x2,y1,y2,s1,s2,N)+R::dnorm(x1_curr,0,sigma,1);
        temp2=LL_id_dep(b1,b2,b3,b4,b01,b02,x1_prev,sigma,mu,theta,x1,x2,y1,y2,s1,s2,N)+R::dnorm(x1_prev,0,sigma,1);
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

double target_theta_dep(double b1, double b2, double b3, double b4, double b01, double b02,
                        NumericVector phi, double mu,NumericVector theta, double rho, double sigma_epsilon,NumericVector x1,NumericVector x2,
                        NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,IntegerVector id,int N,IntegerVector ni, double x1_new, int k,IntegerVector maxni)
{
        double sum = 0.0, func=0.0;
        int m;
        double mu_theta=0.0;
        double var_theta=0.0;

        NumericVector temp=theta;
        for(m=0; m<N; m++)
        {
                if(m==k) {
                        temp[m]=x1_new;
                }
                else
                {
                        temp[m]=theta[m];
                }
        }
        if(maxni[k]>1){
                if(ni[k]>1&&ni[k]<(maxni[k]-1)){
                        mu_theta=(temp[k+1]+temp[k-1])*rho/(1+rho*rho);
                        var_theta=sigma_epsilon/(1+rho*rho);
                }
                else if(ni[k]==1){
                        mu_theta=temp[k+1]*rho;
                        var_theta=sigma_epsilon;
                }
                else {
                        mu_theta=(temp[k-1])*rho;
                        var_theta=sigma_epsilon;
                }
        }
        else {
                mu_theta=0;
                var_theta=sigma_epsilon/(1-rho*rho);
        }

        func=func+(temp[k]-mu_theta)*(temp[k]-mu_theta)/var_theta;
        sum=LL_dep(b1,b2,b3,b4,b01,b02,phi,mu,temp,x1, x2, y1,y2, s1, s2, id,N)-0.5*func;
        return (sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double metroplis_theta_dep(double b1,double b2,double b3, double b4, double b01, double b02,
                           NumericVector phi,double mu,NumericVector theta_curr,
                           double rho, double sigma_epsilon,NumericVector x1,NumericVector x2,NumericVector y1,
                           NumericVector y2,IntegerVector s1,IntegerVector s2,IntegerVector id,int N,IntegerVector ni,
                           double x1_prev,int j,double scale,IntegerVector maxni)

{


        double lgratio=0.0,temp1=0.0, temp2=0.0;
        double u=R::runif(0,1);
        double sd = 1.0;
        double mean = 0.0;
        double rnorm = normalRandom_dep()*sd+mean;
        double x1_curr=x1_prev+scale*rnorm;
        temp1=target_theta_dep(b1,b2,b3,b4,b01,b02,phi,mu,theta_curr, rho, sigma_epsilon, x1,x2,y1,y2,s1,s2,id,N,ni, x1_curr,j,maxni);
        temp2=target_theta_dep(b1,b2,b3,b4,b01,b02,phi,mu,theta_curr,rho, sigma_epsilon, x1,x2,y1,y2,s1,s2,id,N,ni, x1_prev,j,maxni);
        lgratio=temp1-temp2;

        if (lgratio>=log(u))
        {
                theta_curr[j]=x1_curr;
        }
        else
        {
                theta_curr[j]=x1_prev;
        }
        return  theta_curr[j];
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List MH_dep(NumericVector b1,NumericVector b2,NumericVector b3,NumericVector b4,NumericVector b01,
            NumericVector b02,NumericVector phi,NumericVector sigma,NumericVector mu,NumericVector theta,
            NumericVector rho,NumericVector sigma_epsilon,NumericVector x1,NumericVector x2,
            NumericVector y1,NumericVector y2,IntegerVector s1,IntegerVector s2,IntegerVector id,int no_pat_p,
            int n_sample,int N_p,NumericVector scale,IntegerVector ni,IntegerVector maxni){
        int i,j;
        double u;
        double sig_theta,theta_pow,theta_lag,mu_rho,sigma_rho,lower,upper;
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
                temp1=LL_dep(y_new,b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i-1],b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b1[i]=y_new ;
                else b1[i]=b1[i-1];

                u=R::runif(0,1);
                y_new=b2[i-1]+scale[1]*R::rnorm(0,1);
                temp1=LL_dep(b1[i],y_new,b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i],b2[i-1],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b2[i]=y_new ;
                else b2[i]=b2[i-1];

                u=R::runif(0,1);
                y_new=b3[i-1]+scale[2]*R::rnorm(0,1);
                temp1=LL_dep(b1[i],b2[i],y_new,b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i],b2[i],b3[i-1],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b3[i]=y_new ;
                else b3[i]=b3[i-1];

                u=R::runif(0,1);
                y_new=b4[i-1]+scale[3]*R::rnorm(0,1);
                temp1=LL_dep(b1[i],b2[i],b3[i],y_new,b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i],b2[i],b3[i],b4[i-1],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b4[i]=y_new ;
                else b4[i]=b4[i-1];

                u=R::runif(0,1);
                y_new=b01[i-1]+scale[4]*R::rnorm(0,1);
                temp1=LL_dep(b1[i],b2[i],b3[i],b4[i],y_new,b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i],b2[i],b3[i],b4[i],b01[i-1],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b01[i]=y_new ;
                else b01[i]=b01[i-1];


                u=R::runif(0,1);
                y_new=b02[i-1]+scale[5]*R::rnorm(0,1);
                temp1=LL_dep(b1[i],b2[i],b3[i],b4[i],b01[i],y_new,phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                temp2=LL_dep(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i-1],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p);
                lgratio=temp1-temp2;
                if (lgratio>=log(u))  b02[i]=y_new ;
                else b02[i]=b02[i-1];


                u=R::runif(0,1);
                y_new=mu[i-1]+scale[6]*R::rnorm(0,1);
                if(y_new>0){
                        temp1=LL_dep(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,y_new,theta,x1, x2, y1,y2, s1, s2,id,N_p)+R::dgamma(y_new,0.001,1000,1);
                        temp2=LL_dep(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,mu[i-1],theta,x1, x2, y1,y2, s1, s2,id,N_p)+R::dgamma(mu[i-1],0.001,1000,1);
                        lgratio=temp1-temp2;

                        if (lgratio>=log(u))  mu[i]=y_new ;
                        else mu[i]=mu[i-1];
                }
                else mu[i]=mu[i-1];



                for (j = 0; j < no_pat_p; j++)
                {
                        double temp;
                        printf("%d\n",j);
                        x1_prev=phi[j];
                        temp=metroplis_phi_dep(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,sigma[i-1],mu[i],theta,x1,x2,y1,y2,s1,s2,id,no_pat_p,N_p,x1_prev,j,scale[7]);
                        phi[j]=temp;
                        printf("%f\n",phi[j]);
                        phi=center_dep(phi,no_pat_p);
                }



                sig_pat=target_sigma_dep(phi,no_pat_p);
                printf("%f\n",sig_pat);
                sigma[i]=1/R::rgamma(0.001+no_pat_p/2,1/(0.001+0.5*sig_pat));
                printf("%f\n",sigma[i]);
                for (j = 0; j < N_p; j++)
                {
                        double temp;
                        printf("%d\n",j);
                        x1_prev=theta[j];
                        temp=metroplis_theta_dep(b1[i],b2[i],b3[i],b4[i],b01[i],b02[i],phi,mu[i],theta,rho[i-1],sigma_epsilon[i-1],x1,x2,y1,y2,s1,s2,id,N_p,ni,x1_prev,j,scale[8],maxni);
                        theta[j]=temp;
                        printf("%f\n",theta[j]);
                        theta=center_dep(theta,N_p);
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

                theta_pow=theta_pow_sum_dep(theta,N_p,ni);
                theta_lag=theta_lag_sum_dep(theta,N_p,ni);
                mu_rho=theta_lag/theta_pow;
                sigma_rho=sigma_epsilon[i-1]/(theta_pow);


                lower=R::pnorm((-1-mu_rho)/sqrt(sigma_rho),0,1,1,0);
                upper=R::pnorm((1-mu_rho)/sqrt(sigma_rho),0,1,1,0);
                u=R::runif(lower,upper);

                printf("murho is %f\n",mu_rho);
                printf("lower bound is %f\n",lower);
                printf("upper bound is bound is %f\n",upper);
                printf("u is is %f\n",u);

                printf("%d\n",i);
                rho[i]=R::qnorm(u,0,1,1,0)*sqrt(sigma_rho)+mu_rho;

                sig_theta=target_sigma_epsilon_dep(theta,N_p,ni,rho[i]);
                sigma_epsilon[i]=1/R::rgamma(0.001+N_p/2,1/(0.001+0.5*sig_theta));


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
                _["sigma"]= sigma,
                _["sigma_epsilon"]= sigma_epsilon,
                _["mu"]= mu,
                _["rho"]= rho);

}



