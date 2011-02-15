require('gsl') # for hyperg_1F1 in the conditional expectation of the signal

################################################################################

# the canonical gamma MME -- but the cMLE is lower variance and just as fast!
gamma.mme <- function(x) { # {{{
  return(c(shape=(mean(as.matrix(x),na.rm=T)/sd(as.matrix(x),na.rm=T))^2,
           scale=var(as.matrix(x),na.rm=T)/mean(as.matrix(x),na.rm=T)))
} # }}}

## very fast approximation to the full Gamma MLE via Minka (2002) at MS Research
## FIXME: move this to C++
##
gamma.mle <- function(x,w=NULL,niter=100,tol=0.000000001,minx=1) { # {{{

  if( is.null(w) ) w <- rep( 1/length(x), length(x) )
  meanlogx <- weighted.mean(log(pmax(x,minx)), w)
  meanx <- weighted.mean(pmax(x,minx), w)
  logmeanx <- log(meanx)
  a <- a0 <- (0.5/(logmeanx-meanlogx))  # from Minka 2002
  update.a <- function(a) {
    ooa <- 1/a
    1/(ooa+((meanlogx-logmeanx+log(a)-digamma(a))/(((ooa-trigamma(a))*(a**2)))))
  }
  for(i in 1:niter) { # usually converges in under 5 iterations
    a <- update.a(a0)
    if(a0-a < tol) break
    else a0 <- a 
  }
  b <- meanx/a
  # cat('Gamma MLE converged in',i,'iterations\n')
  return(c(shape=a, scale=b))

} # }}}

gamma.mode <- function(par) { # {{{
  ifelse(par['shape'] >= 1, (par['shape']-1)*par['scale'], 0)
} # }}}

gamma.bg <- function(object, channel=NULL, channels=c('Cy3','Cy5')) { # {{{
  if(is.null(channel)) {
    perchannel <- lapply(channels, function(x) gamma.bg(object,x))
    names(perchannel) <- channels
    return(perchannel)
  }
  return(apply(negctls(object, channel), 2, function(z) gamma.mle(z)))
} # }}}

gamma.bg.ebayes <- function(object, channel=NULL, channels=c('Cy3','Cy5')){ #{{{
  
  # pool the gamma parameter distributions to estimate hyperparameters 
  
  # integrate(marginalScale, 0, Inf)
  # integrate(marginalShape, 0, Inf)

} # }}}

## FIXME: don't forget to bgcorrect the non-negative control probes too!!!
##
allelic <- function(x,channel=NULL,allele=NULL,mixture=F,hard=F,parallel=F){#{{{

  ## FIXME: turn this into a defmacro already!!1
  if(!is.null(channel) && tolower(channel) %in% c('cy3','cy5')) { # {{{
    getchan <- match.fun(tolower(channel))
    probes <- getchan(x)
    # }}} 
  } else { # {{{
    probes <- 1:dim(x)[1]
  } # }}}

  # if we fit the fuzzy mixture, add a pi1 attribute
  pi0 <- methylated(x)[probes,] / total.intensity(x)[probes,] # i.e., Beta

  ## FIXME: should we weight pi0 by detection here?
  if(hard) pi0 <- round(pi0) # hard == dichotomized 
  if(mixture) pi0 <- beta_unmix(x, parallel=parallel) # post-EM assignments

  intensities <-list(
    signal=((pi0*methylated(x)[probes,]) + ((1-pi0)*unmethylated(x)[probes,])),
     noise=(((1-pi0)*methylated(x)[probes,]) + (pi0*unmethylated(x)[probes,]))
  )
  
  if(is.null(allele)) return(intensities) 
  else return(intensities[[allele]])

} # }}}

gamma.allelic <- gamma.fg <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('signal','noise'), hard=F) { # {{{

  ## FIXME: replace this nonsense with a macro
  if(is.null(channel)) { # {{{
    perchannel <- lapply(channels, function(channel) {
      gamma.allelic(object, channel=channel, allele=allele, hard=hard)
    })
    names(perchannel) <- channels
    return(perchannel)
  } # }}}
  if(is.null(allele)) { # {{{
    perallele <- lapply(alleles, function(allele) {
      gamma.allelic(object, channel=channel, allele=allele, hard=hard)
    })
    names(perallele) <- alleles
    return(perallele)
  } # }}}

  # unless using hard assignments, we'll just set pi0(y) to beta(y)
  apply( allelic(object,channel=channel,allele=allele,hard=hard),
         # nonnegs(object,channel=channel),  # unsure why this causes problems
         2, gamma.mle )

} # }}}

## FIXME: move this to C++
##
gamma.conditional <- function(total, params) { # {{{

  if(length(total) > 1) sapply(total, gamma.conditional, params=params)
  if(total > (params[3]*params[4])+(10*sqrt(params[3])*params[4])){ # mu+sd.bg
    return( total - (params[3]*params[4]) ) # total - mean(bg)
  } else { 
    g <- params[1] # signal shape
    a <- params[2] # signal scale
    d <- params[3] # bg shape
    b <- params[4] # bg scale
    return( 
      integrate( 
        function(x) {
          # print(paste('Computing 1F1(',g,',',g+d,',',total*((1/b)-(1/a)),')'))
          (exp(x*((1/b)-(1/a)))*(total**(1-g-d))*((total-x)**(d-1))*(x**(g-1)))*
          (1/(beta(g,d)*hyperg_1F1(g, g+d, total*((1/b)-(1/a)), strict=F)))*x
          # print(paste('f(',x,')=',num,'/',den,'=',num/den)) # = num*(1/den)*x
        }, 
        0, total
      )$value # else will return a list with value, abs.error, subdivisions, ...
    )
    # integrate(PrSignalGivenTotal, 0, total)
  }

} # }}}

## FIXME: add a qa step for the remapped beta-mixture scheme or don't use it
gamma.signal <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('signal','noise'), hard=F, parallel=F) { # {{{

  if(is.null(channel)) { # {{{
    perchannel <- lapply(channels, function(channel) {
      gamma.signal(object, channel=channel, allele=allele, hard=hard)
    })
    names(perchannel) <- channels
    return(perchannel)
  } # }}}
  if(is.null(allele)) { # {{{
    perallele <- lapply(alleles, function(allele) {
      gamma.signal(object, channel=channel, allele=allele, hard=hard)
    })
    names(perallele) <- alleles
    return(perallele)
  } # }}}

  params <- rbind(gamma.fg(object, channel, allele), gamma.bg(object, channel))
  rownames(params) <- c('fg.shape','fg.scale','bg.shape','bg.scale')
  totals <- rbind(allelic(object, channel=channel, hard=hard)[[allele]], params)
  pars <- (nrow(totals)-nrow(params)+1):(nrow(totals))

  browser()

  ## should parallelize here 
  signal <- sapply(1:dim(totals)[2], function(j) {
    sapply(totals[-pars,j], gamma.conditional, params=totals[pars,j])
  })

  ## Now reconstruct the methylated and unmethylated bg-corrected signals,
  ## using the adjusted betas to push them towards the "appropriate" place.
  ##
  ## FIXME: either parallelize this or rewrite it in C++
  lapply(signal, function(x) gamma.conditional(x, pars))

  ## recover M/U as pi0*signal + (1-pi0)*noise for each allele

} # }}}

## cross-correlation between replicates and such, for testing normalization
spcor <- function(object, reps=NULL, parallel=FALSE, ... ) { # {{{
  combos <- matrix(0,1,2)
  if(is.null(reps)) reps <- 1:dim(object)[2]
  for(i in 1:length(reps)) {
    for(j in i:length(reps)) { 
      if( i != j ) combos <- rbind(combos, c(reps[i], reps[j]))
    }
  }
  x <- methylated(object)/total.intensity(object)
  if( parallel ) {
    require('multicore')
    results <- c(unlist(mclapply(2:nrow(combos), function(m)
      cor(x[,combos[m,1]], x[,combos[m,2]], method="spearman", use="complete")
    )))
  } else {
    results <- c(unlist(lapply(2:nrow(combos), function(m)
      cor(x[,combos[m,1]],x[,combos[m,2]],use="comp")
    )))
  }
  return(results)
} # }}}

spcor.plot <- function(x, ID=NULL, parallel=TRUE) { # {{{
  if(parallel) require(multicore)
  if(is.null(ID)) ID <- factor(pData(x)[['ID']])
  if(is.null(ID)) stop('Cannot work without IDs!') 
  cols <- c('red','green','brown','blue','orange')
  spcs <- mclapply(levels(ID), function(i) spcor(x,which(pData(x)[['ID']]==i)))
  names(spcs) <- levels(ID)
  plot( density(unlist(lapply(spcs, mean, na.rm=T))), lwd=3,
        main='Replicate correlation (black is overall)',
        ylab='Density', xlab='Spearman correlation' )
  levelses <- levels(ID)
  for(l in 1:nlevels(ID)) {
    lines( density(na.omit(spcs[[levelses[l]]])), 
           col=(l%%nlevels(ID))+1, lty=2 )
  }
} # }}}
