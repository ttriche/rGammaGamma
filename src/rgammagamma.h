#ifndef _rgammagamma_h
#define _rgammagamma_h

#include <cmath>    // lgamma, digamma, etc.
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP rgammagamma_binary_EM(SEXP pi0);
  
NumericVector complement( NumericVector x );

NumericVector beta_xform( NumericVector x, NumericVector w );
  
NumericVector beta_wmle( NumericVector x, NumericVector w ); 

double beta_wll( NumericVector par, NumericVector x, NumericVector w );

#endif
