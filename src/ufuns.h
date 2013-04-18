#ifndef _RcppUtils_UFUNS_H
#define _RcppUtils_UFUNS_H

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// utility with positive and negative consumption
vec u_pos(vec c,vec lab,vec h,double alph,double gamm){
	double mgamm = 1-gamm;
	double imgamm = 1/mgamm;
	double g,z,y;
	vec ret(c);
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab[i] ), mgamm);
		z      = y * h[i];
		g      = pow( c[i], mgamm);
		ret[i] = imgamm * g * z;
	}
	return ret;
}

// quadratic approximation function for case of no work for vector
vec u_neg(vec c,double cuto,vec lab,vec h,double alph,double gamm){
	vec diff;
	vec ret(c);
	double grad,hess,y,z;
	double mgamm = 1-gamm;
	double imgamm = 1/mgamm;
	diff = c - cuto;
	for (int i=0;i<c.n_elem;i++){
		y      = pow( exp( alph * lab[i] ), mgamm);
		z      = y * h[i];
		grad   = pow( cuto, -gamm) * z;
		hess   = -(gamm * grad) / cuto;
		ret[i] = imgamm * pow( cuto, 1-gamm) * z  + (grad * diff[i]) + 0.5 * hess * pow(diff[i],2);
	}
	return ret;
}


//' CRRA utility augmented with discrete labor supply and utility from housing
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
mat ufun_discreteL(mat Res, Rcpp::List par, uvec hsize, vec labor){

	BEGIN_RCPP

	uword n = Res.n_rows;
	uword m = Res.n_cols;

    if (labor.size() != n){
		throw std::logic_error("ufun_discreteL:::error. labor vector not equal rows of cons matrix.");
		return R_NilValue;
	}

	if ( hsize.n_elem != Res.n_rows){
		throw std::logic_error("ufun_discreteL:::error. hsize and Res are not conformable");
		return R_NilValue;
	}

	// extract elements from list par and make some paramters
	double theta     = Rcpp::as<double>(par["theta"]);
	double phival    = Rcpp::as<double>(par["phival"]);
	double mu        = Rcpp::as<double>(par["mu"]);
	double gamma     = Rcpp::as<double>(par["gamma"]);
	double cutoff    = Rcpp::as<double>(par["cutoff"]);
	double alpha     = Rcpp::as<double>(par["alpha"]);

	vec phivals;
	phivals << 0 << phival << 1 << endr;
	vec phivec(hsize.size());
	for (int i=0; i<hsize.n_elem; i++) {
		phivec(i) = phivals( hsize( i ) );
	}

	// initiate return objects as zero
	mat util(Res);
	util.zeros();

	// prepare additive housing premium
	mat phimat = repmat(phivec,1,m);
	mat hfac = exp( theta * phimat);

	mat labmat = repmat(labor,1,m);

	// split Res according to cases:
	//		- pos resouces
	//		- neg resources
	uvec ipos  = find( Res >= cutoff );
	uvec ineg  = find( Res < cutoff );
	vec posres = Res.elem( ipos );
	vec negres = Res.elem( ineg );
	vec posh   = hfac.elem( ipos );
	vec negh   = hfac.elem( ineg );
	vec posl   = labmat.elem( ipos );
	vec negl   = labmat.elem( ineg );

	// calculate utility in each case
	vec upos = u_pos(posres,posl,posh,alpha,gamma);
	vec uneg = u_neg(negres,cutoff,negl,negh,alpha,gamma);

	// reassemble vectors from pos/neg
	util.elem(ipos) = upos;

	if (!(ineg.is_empty())){
		util.elem(ineg) = uneg;
	}

	// add additive premium if houseing is of right size
	util = util + mu * phimat;

	return util;

	END_RCPP
}


#endif
