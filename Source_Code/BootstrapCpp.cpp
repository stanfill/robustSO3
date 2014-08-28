#include <RcppArmadillo.h>   
//#include "/Users/stanfill/Documents/rotationsC/rotations/inst/include/rotations.h"
#include "C:/Users/Sta36z/Documents/rotationsC/rotations/inst/include/rotations.h"


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]] 

// [[Rcpp::export]]
double kappaHatCpp(arma::mat Rbar){
  
  arma::vec eigval;
  arma::svd(eigval,Rbar);
  double lamBar=0.0;
  //eigval.print("Eigs: ");
  
  for(int i=0; i<3;i++){
    lamBar += eigval(i)/3;
  }
  
  return 2*lamBar/(1-lamBar);
  
}



// [[Rcpp::export]]
arma::mat GenQ4Cpp(NumericVector rs){
  RNGScope scope;
  int n = rs.size();
  arma::mat Qs(n,4);
  NumericVector theta(n), phi(n), costheta(n);
  double sinr2, sthetai;
  
  //Used to generate the axis of rotation
  costheta = runif(n,-1,1);
  phi = runif(n,-M_PI,M_PI);
  
  for(int i=0; i<n; i++){
    theta[i] = acos(costheta[i]);
    sinr2 = sin(rs[i]/2);
    sthetai = sin(theta[i]);
    
    Qs(i,0) = cos(rs[i]/2);
    Qs(i,1) = sthetai*cos(phi[i])*sinr2;
    Qs(i,2) = sthetai*sin(phi[i])*sinr2;
    Qs(i,3) = cos(theta[i])*sinr2;
  }
  
  return Qs;
  
}


// [[Rcpp::export]]
arma::rowvec HnBootCpp(arma::mat Rs, int m, int type){
  
  // if type==1 then intrinsic
  // if type!=1 then extrinsic
  
  int n = Rs.n_rows, i;
  arma::mat Rbarels = mean(Rs);
  arma::mat Rbar(3,3), Qs(n,4);
	double kapHat=0.0;
  arma::rowvec Hni(n), Hn(m);
  NumericVector rs(n);
  
	for(i=0;i<9;i++){
			Rbar[i] = Rbarels[i];
	}
  
  /*Estimate the mean and concentration based on the given rotation matrices
  I don't think where the data are centered matters, just the concentration
  Therefore I'll omit Shat and just center the data as id.SO3 (implicitly)*/
  
  //arma::mat Shat = rotations::projectSO3C(Rbar);
  kapHat = kappaHatCpp(Rbar);
  
  /*Generate angles then quaternions based on fitted distribution
  Test for type before jumping into the bootstrap loop, marginally faster then testing
  within the bootstrap loop*/
  
  if(type==1){
    for(i=0;i<m;i++){
      rs = rotations::rcayleyCpp(n, kapHat);
      Qs = GenQ4Cpp(rs);
      Hni = rotations::HnCppIntrinsic(Qs);
      Hn(i)= max(Hni);
    }
  }else{
    for(i=0;i<m;i++){
      rs = rotations::rcayleyCpp(n, kapHat);
      Qs = GenQ4Cpp(rs);
      Hni = rotations::HnCpp(Qs);
      Hn(i)= max(Hni);
    } 
  }
  return Hn;
  
}

