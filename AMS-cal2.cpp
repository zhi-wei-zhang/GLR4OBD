#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>

#include <exception>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;

double min1(double a, double b){
  double z=0;
  if(a>=b){
    z=b;
  }else{
    z=a;
  }
  
  return(z);
}


double max1(double a, double b){
  double z=0;
  if(a>=b){
    z=a;
  }else{
    z=b;
  }
  
  return(z);
}


arma::vec seq_default(long double from, long double to, long unsigned int length_out){
  return arma::linspace(from, to, length_out);
}

IntegerVector seq_int(int a, int b) {
  IntegerVector vec =  seq(a, b);
  return vec;
}


// [[Rcpp::export]]
arma::vec posteriorHE(arma::vec y,
                      arma::vec n,
                      arma::mat p_sample
){
  int ndose;
  ndose = y.size();
  arma::vec like(ndose);
  arma::vec likeli(p_sample.n_rows);
  like.fill(1);
  int i;
  int j;
  for(i=0;i<ndose;i++){
    likeli.fill(1);
    for(j=0;j<ndose;j++){
      likeli=likeli%pow(p_sample.col(j+(i)*ndose),y[j])%pow(1-p_sample.col(j+(i)*ndose),n[j]-y[j]);
    }
    like[i]=mean(likeli);
  }
  like=like/sum(like)+seq_default(1,ndose,ndose)*0.001;
  return like;
}

// [[Rcpp::export]]
arma::vec posteriorHT(arma::vec y,
                      arma::vec n,
                      arma::mat p_sample
){
  int ndose;
  ndose = y.size();
  arma::vec like(ndose+1);
  arma::vec likeli(p_sample.n_rows);
  like.fill(1);
  int i;
  int j;
  for(i=0;i<=ndose;i++){
    likeli.fill(1);
    for(j=0;j<ndose;j++){
      likeli=likeli%pow(p_sample.col(j+(i)*ndose),y[j])%pow(1-p_sample.col(j+(i)*ndose),n[j]-y[j]);
    }
    like[i]=mean(likeli);
  }
  like=like/sum(like);
  return like;
}


// [[Rcpp::export]]
List post_cal2(arma::mat likeliT,
              arma::mat likeliE,
              arma::mat likeliC,
              arma::mat lfT,
              arma::mat lfE,
              arma::mat lfC,
              int ndose){
  int i;
  arma::vec posT(ndose+1);
  arma::vec posE(ndose);
  arma::vec posC(ndose);
  
  likeliT = likeliT%lfT;
  likeliE = likeliE%lfE;
  likeliC = likeliC%lfC;
  
  for(i=0;i<ndose;i++){
    posT[i] = mean(likeliT.col(i));
    posE[i] = mean(likeliE.col(i));
    posC[i] = mean(likeliC.col(i));
  }
  
  posT[ndose] = mean(likeliT.col(ndose));
  posT=posT/sum(posT);
  posE=posE/sum(posE);
  posC=posC/sum(posC);
  
  List result;
  result["likeliT"] = likeliT;
  result["likeliE"] = likeliE;
  result["likeliC"] = likeliC;
  result["posT"] = posT.t();
  result["posE"] = posE.t();
  result["posC"] = posC.t();
  return result;
}