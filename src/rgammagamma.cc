#include "rgammagamma.h"
#define TOLERANCE 0.001
#define MAXITER 25

SEXP rgammagamma_binary_EM( SEXP pi0 ) { // {{{
 
  double ll0, ll1;
  NumericVector x(pi0) ;
  int n = x.size();
  NumericVector piM(n) ;  // intermediate container
  NumericVector piU(n) ;  // intermediate container
  NumericVector pi1(n) ;  // holds results from E-M
  NumericVector abM(2) ;  // parameters for M distn
  NumericVector abU(2) ;  // parameters for U distn
  abM = beta_wmle(x, x); 
  abU = beta_wmle(x, complement(x));
  ll0 = beta_wll(abM, x, x) + beta_wll(abU, x, complement(x));

  for( int i=0; i < MAXITER; i++ ) {

    // expectation of pi1 given abM, abU
    piM = dbeta(x, abM[0], abM[1]);
    piU = dbeta(x, abU[0], abU[1]);
    pi1 = piM / (piM + piU); // vectorized

    // maximum likelihood estimates of abM, abU given pi1
    abM = beta_wmle(x, pi1);
    abU = beta_wmle(x, complement(pi1));

    // compute the loglikelihood given the updated estimates
    ll1 = beta_wll(abM, x, pi1) + beta_wll(abU, x, complement(pi1));

    // break if loglikelihood stops changing
    if( abs(ll1 - ll0) < TOLERANCE ) break;

  }

  return(pi1);

} // }}}

NumericVector complement( NumericVector x ) { // {{{
  
  int n = x.size();
  NumericVector xx(n);
  for( int i=0; i < n; i++ ) xx[i] = 1 - x[i];
  return(xx);

} // }}}

NumericVector beta_xform( NumericVector x, NumericVector w ) { //{{{

  int n = x.size();
  NumericVector x1(n);
  // double s = 0.5;
  double s = weighted_mean(x, w);
  x1 = ((x*(n-1))+s)/n; // vectorized by sugar 
  return x1;

} // }}}

double beta_wll( NumericVector par, NumericVector x, NumericVector w ) { // {{{

  int n = x.size();
  NumericVector zz(n);
  //
  // == (n*lgamma(a+b)) - (n*lgamma(a)) - (n*lgamma(b)) +
  //    ((a-1)*sum(w*log(x))) + ((b-1)*sum(w*log(complement(x))))
  //
  // or just use Rcpp::sugar to roll through dbeta() 
  zz = (w * dbeta( x, par[0], par[1], 1 ) );  // log
  return(sum(zz));

} // }}}

double weighted_mean( NumericVector x, NumericVector w ) { //{{{

  double y = 0.0;
  int n = w.size();
  for( int i=0; i < n; i++ ) {
    y += (x[i]*w[i]);
  }
  return( y / sum(w) );

} // }}}

double weighted_var( NumericVector x, NumericVector w ) { //{{{

  double s2 = 0.0 ;
  int n = w.size() ;
  double xb = weighted_mean(x, w) ;
  for( int i=0; i < n; i++ ) {
    s2 += ((w[i] * pow(x[i]-xb, 2) ) / sum(w)) ;
  }
  return( s2 ) ;

} // }}}

NumericVector beta_wmme( NumericVector x, NumericVector w ) { // {{{

  NumericVector pars(2);
  double xb = weighted_mean(x, w);
  double s2 = weighted_var(x, w);
  pars[1] = xb * (((xb * (1-xb)) / s2) - 1);
  pars[2] = (1-xb) * (((xb * (1-xb)) / s2) - 1);
  return(pars);

} // }}}

NumericVector beta_wmle( NumericVector x, NumericVector w ) { // {{{

  NumericVector par0(2);
  NumericVector par1(2);
  par1 = par0 = beta_wmme( x, w ); 

  // FIXME: Newton-Raphson or L-BFGS-B? 
  //
  // wll( theta ) = (n*lgamma(a+b)) - (n*lgamma(a)) - (n*lgamma(b)) +
  //                ((a-1)*sum(w*log(x))) + ((b-1)*sum(w*log(complement(x))))
  //
  // grad( theta ) = NumericVector(2) = (g1, g2)
  //   g1( theta ) = digamma(a+b) - digamma(a) + (sum(w*log(x))/sum(w))
  //   g2( theta ) = digamma(a+b) - digamma(b) + (sum(w*log(1-x))/sum(w))
  //
  // hess( theta ) = NumericMatrix(2, 2) = ((dg1a, dxab), (dxab, dg2b))
  // dg1a( theta ) = trigamma(a) - trigamma(a + b)
  // dxab( theta ) = -1 * trigamma(a + b)
  // dg1a( theta ) = trigamma(b) - trigamma(a + b)
  //
  // for( int i = 0; i < MAXITER; i++ ) { 
  //   par1 = par0 - (solve(hess(par0))*grad(par0));
  //   if( abs( beta_wll(par1,x,w) - beta_wll(theta0,x,w)) < TOLERANCE ) break;
  // }

  return(par1);

} // }}}
