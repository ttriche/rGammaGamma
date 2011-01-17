beta.mme <- function(x, w=NULL) { # {{{
  if(is.null(w)) w <- rep(1/length(x), length(x))
  xb <- weighted.mean(x, w, na.rm=T)
  s2 <- sum(na.omit(w*((x-xb)**2))) / sum(w[!is.na(x)])
  a <- xb * (((xb*(1-xb))/s2)-1)
  b <- (1-xb) * (((xb*(1-xb))/s2)-1)
  return(c(a=max(a,0), b=max(b,0)))
} # }}}

beta.mode <- function(x, w=NULL) { # {{{
  pars <- beta.mme(x, w)
  if(!all(pars > 1)) return(NaN)
  else return((pars[1] - 1)/(sum(pars)-2))
} # }}}

betaff.mme <- function(x, w=NULL) { # {{{
  pars <- beta.mme(x, w)
  return(c(mu=pars[1]/(sum(pars)), phi=sum(pars)))
} # }}}

beta.unmix <- function(x, niter=25, parallel=FALSE, tol=0.000001) { # {{{ 

  require(impute)
  if(is(x, 'MethyLumiSet') || is(x, 'MethyLumiM')) x <- betas(x)
  if(anyMissing(x)) x <- impute.knn(x)$data
  if(is(x, 'matrix')) { # {{{
    if(parallel) {
      require(multicore) ## FIXME: OpenMP would be better
      data.matrix(as.data.frame(mclapply(1:dim(x)[2], function(i) {
        beta.unmix(x[,i])
      })))
    } else {
      apply(x, 2, beta.unmix)
    }
  } # }}}

  pi0 <- pi1 <- x
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
      return(list(pi1=pi1, theta0=theta0, theta1=theta1))
    } else { 
      negloglik0 <- negloglik1
      theta0 <- theta1 
      pi0 <- pi1
    }
  }
  if(i == niter) warning("beta.unmix: Failure to converge!")
  return( pi1 )

} # }}}

