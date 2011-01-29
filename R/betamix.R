## faster than the MLE and typically lower variance as well; accepts weights
##
beta.mme <- function(x, w=NULL, ...) { # {{{
  if(any(is.na(x))) stop("Cannot handle NA values!")
  if(is.null(w)) w <- rep(1/length(x), length(x))
  xb <- weighted.mean(x, w, na.rm=T)
  s2 <- sum(w*((x-xb)**2)) / sum(w)
  a <- xb * (((xb*(1-xb))/s2)-1)
  b <- (1-xb) * (((xb*(1-xb))/s2)-1)
  mme <- c(a=max(a,0), b=max(b,0))
  return(mme)
} # }}}

## from Smithson & Verkuilen 2006; shrink towards either 0.5 or the mean
##
beta.transform <- function(x, w=NULL, to.mean=FALSE, to.mode=TRUE, s=0.5){ #{{{

  stopifnot( (max(x)<=1) && (min(x)>=0) )
  n <- length(x)
  if( is.null(w) ) w <- rep(1, n)
  else w <- (w / (sum(w)/n))
  if(to.mean) s <- weighted.mean(x, w, na.rm=T)
  if(to.mode) s <- pmax(0.01, pmin(0.99, beta.mode(beta.mme(x, w))))
  return( ((x*(n-1))+s) / n ) 

} # }}}

## this should be farmed out to C++ (actually the entire unmixing should be)
##
beta.mle <- function(x, w=NULL, transform=TRUE) { #{{{

  n <- length(x)
  if(is.null(w)) w <- rep(1, n)
  w <- w / (sum(w)/n)
  par0 <- beta.mme(x, w) 
  if(transform) x <- beta.transform(x, w)
  fn <- function(params) -1 * sum(w * dbeta(x, params[1], params[2], log=T))
  gr <- function(params) {  # {{{
    a <- params[1]
    b <- params[2]
    return( c( a=(digamma(a+b) - digamma(a) + (sum(w*log(x))/n)),
               b=(digamma(a+b) - digamma(b) + (sum(w*log(1-x))/n)) ) )
  } # }}}
  res <- try(optim(par0, fn, gr, method="L-BFGS-B", lower=0.1, upper=2**16))
  if(class(res) == 'try-error') {
    return(par0)
  } else {
    return(res$par)
  }

} # }}}

## can be used to set the priors for the mixture model a bit more effectively
## 
beta.mode <- function(a, b) { # {{{
  if(length(a) > 1) {
    b <- a[2]
    a <- a[1]
  }
  if(all(c(a,b) < 1)) return(NaN)
  if(a < 1) return(1)
  if(b < 1) return(0)
  return( (a-1) / (a+b-2) )
} # }}}

## beta mixture model for a "more empirical" prior on signal/noise guesses
## farm this out to C++ to avoid irritation from R being slow!
## 
beta.unmix <- function(x, niter=50, tol=0.01) { # {{{ 

  require(multicore) # until I move this to C++
  if(is(x, 'MethyLumiSet') || is(x, 'MethyLumiM')) {
    x <- pmax(1, methylated(x)) / pmax(2, (methylated(x)+unmethylated(x)))
    return(apply(x, 2, beta.unmix))
  }
  if(anyMissing(x)) stop("NA values are not supported; fix with impute.knn")

  pi1 <- pi0 <- list( M=x, U=1-x )
  theta <- mclapply( pi0, function(pi0) beta.mme(x, pi0) )
  ## theta <- mclapply( pi0, function(pi0) beta.mle(x, pi0) )
  getll <- function(theta, pi0, x) { # {{{
    sum(pi0 * pmin(10000, dbeta(x, theta[1], theta[2], log=T)))
  } # }}}
  sumll <- function(theta, pi0, x) { # {{{
    sum(unlist(lapply(names(pi0), function(nm) {
      getll(theta[[nm]], pi0[[nm]], beta.transform(x, pi0[[nm]]))
    })))
  } # }}}
  ll0 <- sumll( theta, pi0, x )

  ## FIXME: do this in C++
  for(i in 1:niter) {
    if(options()$verbose) cat("Iteration", i, "...\n")
    pi1 <- mclapply(theta, function(th) dbeta(x, th[1], th[2]))
    pi0 <- mclapply(pi1, function(wts) pmax(0.01,wts)/pmax(0.02,(pi1$M+pi1$U)))
    theta <- mclapply(pi0, function(wts) beta.mme(x, wts))
    ll <- sumll(theta, pi0, x)
    discrep <- ll - ll0
    if( abs(discrep) < tol ) return( pi0$M )
    else ll0 <- ll
  }
  warning("Failed convergence: discrepancy ",discrep," at tolerance ",tol)
  return( pi0$M ) 

} # }}}

