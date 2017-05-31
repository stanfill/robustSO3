#include <RcppArmadillo.h>   
#include "/Users/stan070/Documents/Rotations/rotationsC/rotations/inst/include/rotations.h"
//#include "/Users/stanfill/Documents/rotationsC/rotations/inst/include/rotations.h"
//#include "C:/Users/Sta36z/Documents/rotationsC/rotations/inst/include/rotations.h"


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]
arma::mat QsBootstrap(arma::mat Qs){
  RNGScope scope;
  
  int n = Qs.n_rows;
  int m = Qs.n_cols;
  arma::mat QsBoot(n,m);
  NumericVector rows = runif(n,0,n);
  rows = floor(rows);

  for(int i=0;i<n;i++){
    QsBoot.row(i) = Qs.row(rows(i));
  }
  
  return QsBoot;
}

// [[Rcpp::export]]
arma::rowvec HnNonParaBootCpp(arma::mat Qs, int m, int type){
  /*rangle is the function used to generate the misorientation angles*/
  
  // if type==1 then intrinsic
  // if type!=1 then extrinsic
  
  int n = Qs.n_rows,i;
  arma::rowvec Hni(n), Hn(m);
  arma::mat QsResamp(n,4);
  
  /*Generate angles then quaternions based on fitted distribution
  Test for type before jumping into the bootstrap loop, marginally faster then testing
  within the bootstrap loop*/
  
  if(type==1){
    for(i=0;i<m;i++){

      Hni = rotations::HnCppIntrinsic(QsBootstrap(Qs));
      Hn(i)= max(Hni);
    }
  }else{
    for(i=0;i<m;i++){

      Hni = rotations::HnCpp(QsBootstrap(Qs));
      Hn(i)= max(Hni);
    } 
  }
  return Hn;
  
}

