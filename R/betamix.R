## this is the precision term for beta regression
## 
beta.phi <- function(x, w=NULL) { # {{{
  return(sum(beta.mme(x, w)))
} # }}}

## faster than the MLE and typically lower variance as well; accepts weights
##
beta.mme <- function(x, w=NULL) { # {{{
  if(is.null(w)) w <- rep(1/length(x), length(x))
  xb <- weighted.mean(x, w, na.rm=T)
  s2 <- sum(na.omit(w*((x-xb)**2))) / sum(w[!is.na(x)])
  a <- xb * (((xb*(1-xb))/s2)-1)
  b <- (1-xb) * (((xb*(1-xb))/s2)-1)
  return(c(a=max(a,0), b=max(b,0)))
} # }}}

## used to set the priors for the mixture model a bit more effectively
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
## 
beta.unmix <- function(x, niter=100, rows=F, parallel=F, tol=0.01) { # {{{ 

  if(is(x, 'MethyLumiSet') || is(x, 'MethyLumiM')) { # {{{
    x <- pmax(1, methylated(x)) / pmax(2, (methylated(x)+unmethylated(x)))
    x <- pmax(0.000001, pmin(0.999999, x)) # squeeze play 
  } # }}}
  if(anyMissing(x)) stop("NA values are not supported; fix with impute.knn")
  if(is(x, 'matrix')) { # {{{
    if(parallel) {
      require(multicore) ## FIXME: OpenMP would be better
      if(rows) {
        t(data.matrix(as.data.frame(mclapply(1:dim(x)[1], function(i) {
          beta.unmix(x[i,])
        }))))
      } else {
        data.matrix(as.data.frame(mclapply(1:dim(x)[2], function(i) {
          beta.unmix(x[i,])
        })))
      }
    } else {
      if(rows) apply(x, 1, beta.unmix)
      else apply(x, 2, beta.unmix)
    }
  } # }}}

  pi0 <- pmax(0.01, pmin(0.99, x))
  theta0 <- list(U=beta.mme(x, 1-pi0), M=beta.mme(x, pi0))

  ## FIXME: do in C++
  negloglik0 <- 10000
  for(i in 1:niter) {
    pi1 <- dbeta(x, theta0$M[1], theta0$M[2])
    pi1 <- pi1 / (pi1 + dbeta(x, theta0$U[1], theta0$U[2]))
    theta1 <- list(U=beta.mme(x, 1-pi1), M=beta.mme(x, pi1))
    negloglik1 <- -1 * sum(log(((1-pi1)*dbeta(x, theta1$U[1], theta1$U[2])) +
                               (pi1*dbeta(x, theta1$M[1], theta1$M[2])))) 
    if( abs(negloglik1 - negloglik0) < tol ) {
      return( pi1 )
    } else { 
      discrep <- negloglik1-negloglik0
      negloglik0 <- negloglik1
      theta0 <- theta1 
      pi0 <- pi1
    }
  }
  if(i == niter) {
    warning("Failed convergence: discrepancy ", discrep, " at tolerance ",tol)
  }
  return( pi1 )

} # }}}

