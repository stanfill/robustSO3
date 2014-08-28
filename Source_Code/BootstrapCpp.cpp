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

//Necessary but uninteresting functions
arma::rowvec RdistCArma(arma::mat Q1, arma::rowvec Q2){
  /*Compute the geodesic distance between quaternions Q1 and Q2*/
  /* Q1 must be an n-by-4 matrix with quaternion rows and Q2 a single (1x4) quaternion*/
  
  int n = Q1.n_rows, i=0; 
	double cp;
	arma::rowvec rs(n);
	
	for(i=0;i<n;i++){
		
		cp = sum(Q1.row(i)*Q2.t());
		rs(i) = acos(2*cp*cp-1);
		
	}
	
	return rs;
}

arma::rowvec HnCpp(arma::mat Qs){
  //Compute the Hn tests statistics
  
  int n = Qs.n_rows, i=0;
  arma::mat T = Qs.t()*Qs;
  arma::mat eigvec, eigvecJ;
  arma::vec eigval, eigvalJ;
  arma::eig_sym(eigval,eigvec,T);
  arma::rowvec Hn(n);
  arma::rowvec Qj;
  arma::mat Tj;

  for(i = 0;i<n; i++){
    Qj = Qs.row(i);
    
    Tj = T-Qj.t()*Qj;
    arma::eig_sym(eigvalJ,eigvecJ,Tj);
    Hn(i)=(n-2)*(1+eigvalJ(3)-eigval(3))/(n-1-eigvalJ(3));
    
  }
  return Hn;
}

// [[Rcpp::export]]
arma::rowvec HnCppIntrinsic(arma::mat Qs){
  
  //Compute the intrinsic Hn tests statistics
  
  int n = Qs.n_rows, i=0, j=0;

  //Get T matrix of whole sample to make it easier later on
  arma::mat T = Qs.t()*Qs;
  arma::rowvec Qhat = rotations::meanQ4C(Qs);
  arma::rowvec dists(n);

  //Sum of squared geometric distances between proj. mean and each obs
  double SSE = 0.0, SSEJ=0.0;
  dists = square(RdistCArma(Qs,Qhat));
  SSE = sum(dists);

  arma::rowvec Hn(n);
  arma::rowvec Qhatj;
  arma::mat QsJ(n,4);
  
  //Variables for reduced sample mean
  arma::rowvec Qj, distsJ(n-1);
  arma::mat Tj(4,4), eigvecJ(4,4);

  for(i = 0;i<n; i++){
    
    QsJ.resize(n,4);
    
    for(j = 0; j<n; j++){
      if(j!=i){
        QsJ.row(j) = Qs.row(j);
      }
    }
    
    QsJ.shed_row(i);
    
    //Compute projected mean when jth row is cut out
    Qhatj = rotations::meanQ4C(QsJ);
    distsJ = square(RdistCArma(QsJ,Qhatj));
    SSEJ = sum(distsJ);
    
    //Rcpp::Rcout << "SSEJ: " << SSEJ << std::endl;
    
    Hn(i)=(n-2)*(SSE-SSEJ)/(SSEJ);
    
  }

  return Hn;
  
}

NumericVector rcayleyCpp(int n, double kappa){
  RNGScope scope;
  NumericVector bet(n), alp(n), theta(n);
  
  bet = rbeta(n,kappa+0.5,1.5);
  alp = rbinom(n,1,0.5);
  
  for(int i = 0; i<n; i++){
    
    theta[i] = acos(2*bet[i]-1)*(1-2*alp[i]);
  
  }
  return theta;
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
      rs = rcayleyCpp(n, kapHat);
      Qs = GenQ4Cpp(rs);
      Hni = HnCppIntrinsic(Qs);
      Hn(i)= max(Hni);
    }
  }else{
    for(i=0;i<m;i++){
      rs = rcayleyCpp(n, kapHat);
      Qs = GenQ4Cpp(rs);
      Hni = HnCpp(Qs);
      Hn(i)= max(Hni);
    } 
  }
  return Hn;
  
}

