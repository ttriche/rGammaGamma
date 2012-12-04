require('gsl') # for hyperg_1F1 in the conditional expectation of the signal

################################################################################

## the canonical gamma MME -- but the cMLE is lower variance and just as fast!
gamma.mme <- function(x) { # {{{
  return(c(shape=(mean(as.matrix(x),na.rm=T)/sd(as.matrix(x),na.rm=T))^2,
           scale=var(as.matrix(x),na.rm=T)/mean(as.matrix(x),na.rm=T)))
} # }}}

## fast approximation to the full Gamma MLE via Minka (2002) at MS Research
gamma.mle <- function(x,w=NULL,niter=100,tol=0.1,minx=1) { # {{{

  if( is.null(w) ) w <- rep( 1, length(x) )
  meanlogx <- weighted.mean(log1p(x), w)
  meanx <- weighted.mean(x, w)
  logmeanx <- log(meanx)
  a <- a0 <- (0.5/(logmeanx-meanlogx))  # from Minka 2002
  if(is.nan(a)) stop('NaN starting estimate')
  update.a <- function(a) {
    ooa <- 1/a
    1/(ooa+((meanlogx-logmeanx+log(a)-digamma(a))/(((ooa-trigamma(a))*(a**2)))))
  }
  for(i in 1:niter) { # usually converges in under 5 iterations
    a <- update.a(a0)

    if(abs(a0-a) < tol) break
    else a0 <- a 
  }
  b <- meanx/a
  # cat('Gamma MLE converged in',i,'iterations\n')
  return(c(shape=a, scale=b))

} # }}}

## FIXME: move this to C++ as soon as humanly possible (ideally with matrix arg)
gamma.integral <- function(total, params, offset=50, minx=1) { # {{{

  ## this bit is the most obvious "farm me out to C++" piece of all...
  if(length(total) > 1) return(sapply(total, gamma.integral, params=params))
  
  g = params[1]
  a = params[2]
  d = params[3]
  b = params[4]
  bg.mean = d * b 
  bg.sd = sqrt( d * b * b )
  ## cat('total =',total,'... bg.mean =',bg.mean,'... bg.sd =',bg.sd,"\n")
  if(total > ( bg.mean + ( 3 * bg.sd ) )) {
    return(pmax(total - bg.mean, minx))
  } else {
    ## FIXME: need to write this as a function object for C++ to integrate it
    res = try(
      integrate( 
        function(x) {
          (exp(x*((1/b)-(1/a)))*(total**(1-g-d))*((total-x)**(d-1))*(x**(g-1)))*
          (1/(beta(g,d)*hyperg_1F1(g, g+d, total*((1/b)-(1/a)), strict=F)))*x 
        }, 
        0, total
      )$value # else will return a list with value, abs.error, subdivisions, ...
    ) # i.e., integrate(PrSignalGivenTotal, /* from */ 0, /* to */ total);
    if(class(res) == 'try-error') {
      return(pmax(total-bg.mean, minx)+offset)
    } else {
      return(pmax(res, minx)+offset)
    }
  }

} # }}}
