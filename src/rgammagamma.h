#ifndef _rgammagamma_h
#define _rgammagamma_h

#include <cmath>
#include <Rcpp.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>

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

extern "C" double beta( double a, double b );

extern "C" double digam( double x );

extern "C" double trigam( double x );

extern "C" double hyperg1f1( double a, double b, double x );

NumericVector gamma_mle( NumericVector x );

NumericVector gamma_wmle( NumericVector x, NumericVector w );

double gamma_conv( double x, NumericVector params );

double gslfunc_example(double x, void * params);

double gslint_example(double x, double t, void * params);

double convfn(double x, double t, double g, double a, double d, double b);

NumericMatrix gamma_conv( NumericMatrix x, NumericMatrix params );

#endif
