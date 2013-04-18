

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "ufuns.h"

using namespace Rcpp

//' C++ module for computation of V = U + b*EV
//'
//' same as other util functions but with discrete labor supply choice. no time separability, i.e. labor supply is not implied.
//' @param Resources numeric matrix of consumption levels
//' @param hsize numeric vector of house sizes
//' @param labor numeric vector of discrete labor supply choices
//' @param params list of parameters 
//' \itemize{
//' \item{theta}{elasticity of substitution between c and h}
//' \item{phival}{value of relative utility difference flat vs house}
//' \item{mu}{weight on additive utility premium}
//' \item{gamma}{coefficient of relative risk aversion}
//' \item{cutoff}{minimum level of consumption. below cutoff, u(c) is quadratically approximated.}
//' \item{alpha}{coefficient on labor}
//' }
//' @return numeric matrix of utility values 
// [[Rcpp::Export]]
List util_module(NumericMatrix cashR, NumericMatrix saveR, NumericMatrix EVR, NumericVector hsizeR, NumericVector laborR, List par){

	int n = cashR.nrow();   // number of states
	int m = cashR.ncol();	// number of discrete choices
	int k = saveR.ncol();   // number of savings choices

	// some input checks
	if (n != saveR.nrow()){
		throw std::logic_error( "util_module: saveR not equal rows as cashR");
		return R_NilValue;
	} 
	if (n != EVR.nrow()){
		throw std::logic_error( "util_module: EVR not equal rows as cashR");
		return R_NilValue;
	} 
	if (k != EVR.ncol()){
		throw std::logic_error( "util_module: EVR not equal cols as saveR");
		return R_NilValue;
	} 

	// map to R to arma
	arma::mat cash(cashR.begin(), n, m, false);	// advanced constructor. no copying.
	arma::mat EV(EVR.begin(), n, k, false);	
	arma::mat save(saveR.begin(), n, k, false);	
	arma::vec labor(laborR.begin(),n,false);
	arma::uvec hsize(hsizeR.begin(),n,false);
	// Rcpp::List par( par_ ) ;

	// allocate tmpcash, utility, and W matrices
	arma::fmat::fixed<n,k> tmpcash;
	arma::fmat::fixed<n,k> util;
	arma::fmat::fixed<n,k> W;

	// allocate out matrices
	arma::fmat::fixed<n,m> rety;
	arma::umat::fixed<n,m> retiy;

	// allocate for maximization loop
	arma::uword j;
	arma::uvec::fixed<n> iy;
	arma::vec::fixed<n> y;
	arma::rowvec::fixed<k> tmpvec;

	
	// loop over discrete choices
	for (int i=0;i<m;i++){

		// zero out iy and y
		iy.zeros();
		y.zeros();

		// fill tmp with copies of cash
		tmpcash = repmat(cash.col(i),1,k);
		// get consumption at each savings choice
		tmpcash = tmpcash - save;
		// get utility
		util    = ufun_discreteL( tmpcash, par, hsize, labor );
		// get lifetime utility
		W       = util + EV;

		// maximize
		// loop over rows of W and find maximal value and it's index for each row.
		for (int i=0; i<n; i++) {
			tmpvec = W.row(i);
			y(i) = tmpvec.max( j );
			iy(i) = j + 1; // go to 1 based indices
		}
		// put into return objects
		rety.col(i)  = y;
		retiy.col(i) = iy;
	}

	Rcpp::List list = Rcpp::List::create( _["values"] = rety, _["indices"] = retiy );
	return list;
}


