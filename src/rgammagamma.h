#ifndef _rgammagamma_h
#define _rgammagamma_h

#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_sf.h>
#include "DEIntegrator.h"

using namespace Rcpp;
using namespace std;

RcppExport SEXP rgammagamma_binary_EM(SEXP pi0);

RcppExport SEXP rgammagamma_ternary_EM(SEXP pi0);

RcppExport SEXP rgammagamma_gamma_mle( SEXP x );

RcppExport SEXP rgammagamma_gamma_wmle( SEXP x, SEXP w );

NumericVector complement( NumericVector x );

NumericVector beta_xform( NumericVector x, NumericVector w );
  
NumericVector beta_wmle( NumericVector x, NumericVector w ); 

NumericVector beta_wmme( NumericVector x, NumericVector w ); 

double beta_wll( NumericVector par, NumericVector x, NumericVector w );

double weighted_mean( NumericVector x, NumericVector w );

double weighted_var( NumericVector x, NumericVector w );

extern "C" double digam( double x );

extern "C" double trigam( double x );

NumericVector gamma_mle( NumericVector x );

NumericVector gamma_wmle( NumericVector x, NumericVector w );

double gamma_conv( double x, NumericVector params );

NumericMatrix gamma_conv( NumericMatrix x, NumericMatrix params );

// problematic due to type casting
/* class intfn { 
  public:
  double operator()(double t, double g, double a, double d, double b) const {
    return(
    // note the GSL fucntions for hyperg_1f1 and beta below
    //
    //(exp(x*((1/b)-(1/a)))*(t**(1-g-d))*((t-x)**(d-1))*(x**(g-1)))*
    //(1/(beta(g,d)*hyperg_1F1(g, g+d, t*((1/b)-(1/a)), strict=F)))*x
    );
  }
}; */

#endif
