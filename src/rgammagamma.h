#ifndef _rgammagamma_h
#define _rgammagamma_h

#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_sf.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP rgammagamma_binary_EM(SEXP pi0);

RcppExport SEXP rgammagamma_ternary_EM(SEXP pi0);

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

NumericVector gamma_wmle( NumericVector x, NumericVector w );

#endif
