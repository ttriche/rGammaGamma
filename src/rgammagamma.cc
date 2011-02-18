#include "rgammagamma.h"
#define TOLERANCE 0.0001
#define MAXITER 25

RCPP_MODULE(gamma) { // {{{
  function("mle", &gamma_mle );
  function("wmle", &gamma_wmle );
} // }}}

RCPP_MODULE(beta) { // {{{
  function("xform", &beta_xform);
  function("wmme", &beta_wmme);
  function("wmle", &beta_wmle);
  function("wll", &beta_wll);
} // }}}

RcppExport SEXP rgammagamma_binary_EM( SEXP pi0 ) { // {{{
 
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

RcppExport SEXP rgammagamma_ternary_EM( SEXP pi0 ) { // {{{
 
  double ll0, ll1;
  NumericVector x(pi0) ;
  int n = x.size();
  NumericVector piM(n) ;  // intermediate container
  NumericVector piH(n) ;  // intermediate container
  NumericVector piU(n) ;  // intermediate container
  NumericVector pi1(n) ;  // holds results from E-M
  NumericVector abM(2) ;  // parameters for M distn
                          // the uniform component just gets (a,b) == (1,1)
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

RcppExport SEXP rgammagamma_gamma_wmle( SEXP x, SEXP w ) { // {{{

  NumericVector xx(x);
  NumericVector ww(w);
  return( Rcpp::wrap(gamma_wmle( xx, ww )) );

} // }}}

RcppExport SEXP rgammagamma_gamma_mle( SEXP x ) { // {{{

  NumericVector xx(x);
  return( Rcpp::wrap(gamma_mle(xx)) );

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

extern "C" double digam( double x ) { // {{{
  return(gsl_sf_psi(x));  // this is obnoxious
} // }}}

extern "C" double trigam( double x ) { // {{{
  return(gsl_sf_psi_1(x)); // this is obnoxious
} // }}}

NumericVector gamma_wmle( NumericVector x, NumericVector w ) { // {{{

  // no NAs allowed!
  NumericVector pars(2);

  double mlnx = weighted_mean( log(pmax(x,1)), w );
  double mx = weighted_mean( pmax(x,1), w ); 
  double lnmx = log( mx );
  double a = 0.5/(lnmx-mlnx);
  double lna = log(a); 
  double z = a;

  try {
    // from Minka 2002
    for( int i = 1; i < MAXITER; i++ ) { // usually in under 5 iterations
      z = 1/((1/a)+((mlnx-lnmx+log(a)-digam(a))/(((1/a)-trigam(a))*(a*a))));
      if( abs(z - a) < TOLERANCE ) break; 
      else a = z;
    }
    pars[0] = a;
    pars[1] = mx/a;
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)" );
  }
  return(pars);

} // }}}

NumericVector gamma_mle( NumericVector x ) { // {{{

  // no NAs allowed!
  NumericVector pars(2);
  int n = x.size();

  double mlnx = sum(log(pmax(x,1)))/n;
  double mx = sum(pmax(x,1))/n;
  double lnmx = log(mx);
  double a = 0.5/(lnmx-mlnx);
  double lna = log(a); 
  double z = a;

  try {
    // from Minka 2002
    for( int i = 1; i < MAXITER; i++ ) { // usually in under 5 iterations
      z = 1/((1/a)+((mlnx-lnmx+log(a)-digam(a))/(((1/a)-trigam(a))*(a*a))));
      if( abs(z - a) < TOLERANCE ) break; 
      else a = z;
    }
    pars[0] = a;
    pars[1] = mx/a;
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)" );
  }
  return(pars);

} // }}}

// project: move the conditional expectation integral from R to C++ (here)
double gamma_conv( double x, NumericVector params ) { // {{{

  double z;
  if(x > (params(3)*params(4))+(3*sqrt(params(3))*params(4))) { // mu+sd.bg
    z = x - (params(3)*params(4)); // just subtract the background mean
  } else {
    // intfn fn;
    int evals;
    double err;
    // std::setprecision(10);
    // z = DEIntegrator<intfn>::Integrate(fn, 0, x, 1e-6, evals, err);
  }
  return(z);

} // }}}

// project: stride through the columns/rows of a signal matrix for the integral
NumericMatrix gamma_conv( NumericMatrix x, NumericMatrix params ) { // {{{

  // walk through 'x' and 'params' 
  for( int i = 0; i < x.ncol(); i++ ) {
    for( int j = 0; i < x.nrow(); j++ ) {
      x((j+1), (i+1)) = gamma_conv(x((j+1),(i+1)), params(_, (i+1)));
    }
  }
  return(x);

} // }}}
