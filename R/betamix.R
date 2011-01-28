## this is the precision term for beta regression
## 
beta.phi <- function(x, w=NULL) { # {{{
  return(sum(beta.mme(x, w)))
} # }}}

## faster than the MLE and typically lower variance as well; accepts weights
##
beta.mme <- function(x, w=NULL, ...) { # {{{
  if(is.null(w)) w <- rep(1/length(x), length(x))
  xb <- weighted.mean(x, w, na.rm=T)
  s2 <- sum(na.omit(w*((x-xb)**2))) / sum(w[!is.na(x)])
  a <- xb * (((xb*(1-xb))/s2)-1)
  b <- (1-xb) * (((xb*(1-xb))/s2)-1)
  return(c(a=max(a,0), b=max(b,0)))
} # }}}

## from Smithson & Verkuilen 2006; shrink towards either 0.5 or the mean
##
beta.transform <- function(x, w=NULL, to.mean=TRUE, s=0.5) { # {{{

  stopifnot( (max(x)<=1) && (min(x)>=0) )
  n <- length(x)
  if( is.null(w) ) w <- rep(1, n)
  else w <- (w / (sum(w)/n))
  if(to.mean) s <- weighted.mean(x, w, na.rm=T)
  ((x*(n-1))+s) / n

} # }}}

## this should be farmed out to C++ (actually the entire unmixing should be)
##
beta.mle <- function(x, w=NULL, transform=TRUE, to.mean=TRUE) { #{{{

  n <- length(x)
  if(is.null(w)) w <- rep(1, n)
  w <- w / (sum(w)/n)
  cat(length(w))
  par0 <- beta.mme(x, w) 
  if(transform) x <- beta.transform(x, w)
  cat(length(w))
  fn <- function(params) -1 * sum(w * dbeta(x, params[1], params[2], log=T))
  grad <- function(params) {  # {{{
    a <- params[1]
    b <- params[2]
    c( a=(digamma(a+b) - digamma(a) + (sum(w*log(x))/n)),
       b=(digamma(a+b) - digamma(b) + (sum(w*log(1-x))/n)) )
  } # }}}
  res <- try(optim(par0, fn, grad, method="L-BFGS-B", lower=0.0001))
  if(class(res) == 'try-error') 
    return(par0)
  else 
    return(res$par)

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
beta.unmix <- function(x, niter=25, tol=0.001) { # {{{ 

  if(is(x, 'MethyLumiSet') || is(x, 'MethyLumiM')) {
    x <- pmax(1, methylated(x)) / pmax(2, (methylated(x)+unmethylated(x)))
  }
  if(anyMissing(x)) stop("NA values are not supported; fix with impute.knn")

  pi0.0 <- pi0 <- list( M=beta.transform(x), U=1-beta.transform(x) )
  theta <- lapply( pi0, function(pi0) beta.mle(x, pi0) )
  getll <- function(pi0, theta) sum(pi0 * dbeta(x, theta[1], theta[2], log=T))
  loglik0 <- sum(sapply(names(pi0), function(al) getll(theta[[al]], pi0[[al]])))

  ## FIXME: do this in C++
  for(i in 1:niter) {
    cat("Iteration", i, "...\n")
    pi1 <- lapply(theta, function(params) dbeta(x, params[1], params[2]))
    browser()
    pi0 <- pi1$M / (pi1$M + pi1$U)
    theta <- lapply( names(pi0), function(allele) beta.mle(x, pi0[[allele]]))
    loglik <- sum(sapply(names(pi0),function(al) getll(theta[[al]],pi0[[al]])))
    discrep <- loglik - loglik0  
    if( abs(discrep) < tol ) return( pi0 )
    else loglik0 <- loglik
  }
  warning("Failed convergence: discrepancy ",discrep," at tolerance ",tol)
  return( pi0 )

} # }}}

