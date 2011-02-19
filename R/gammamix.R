require('gsl') # for hyperg_1F1 in the conditional expectation of the signal

################################################################################

## the canonical gamma MME -- but the cMLE is lower variance and just as fast!
gamma.mme <- function(x) { # {{{
  return(c(shape=(mean(as.matrix(x),na.rm=T)/sd(as.matrix(x),na.rm=T))^2,
           scale=var(as.matrix(x),na.rm=T)/mean(as.matrix(x),na.rm=T)))
} # }}}

## fast approximation to the full Gamma MLE via Minka (2002) at MS Research
gamma.mle <- function(x,w=NULL,niter=100,tol=0.000000001,minx=1) { # {{{

  if( is.null(w) ) w <- rep( 1, length(x) )
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
    if(abs(a0-a) < tol) break
    else a0 <- a 
  }
  b <- meanx/a
  # cat('Gamma MLE converged in',i,'iterations\n')
  return(c(shape=a, scale=b))

} # }}}

## faster approximation (in C++)
gamma.cmle <- function(x,w=NULL) { # {{{

  if( !is.null(w) ) .Call('rgammagamma_gamma_wmle', x, w)
  else .Call('rgammagamma_gamma_mle', x)

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
allelic <- function(x,channel=NULL,allele=NULL,mixture=F,hard=F,parallel=F){#{{{

  ## FIXME: turn this into a defmacro already!!1
  if(!is.null(channel)) {
    probes <- getProbesByChannel(x, channel)
  } else { 
    probes <- 1:dim(x)[1]
  } 

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

  apply( allelic(object,channel=channel,allele=allele), 2, gamma.mle )
         # nonnegs(object,channel=channel),  # unsure why this causes problems

} # }}}

## FIXME: move this to C++ as soon as humanly possible (ideally with matrix arg)
gamma.conditional <- function(total, params, minx=1) { # {{{

  if(length(total) > 1) sapply(total, gamma.conditional, params=params)
  if(total > (params[3]*params[4])+(3*sqrt(params[3])*params[4])){ # mu+sd.bg
    return( total - (params[3]*params[4]) ) # total - mean(bg)
  } else { 
    g <- params[1] # signal shape
    a <- params[2] # signal scale
    d <- params[3] # bg shape
    b <- params[4] # bg scale
    res <- try( 
      integrate( 
        function(x) {
          # print(paste('Computing 1F1(',g,',',g+d,',',total*((1/b)-(1/a)),')'))
          (exp(x*((1/b)-(1/a)))*(total**(1-g-d))*((total-x)**(d-1))*(x**(g-1)))*
          (1/(beta(g,d)*hyperg_1F1(g, g+d, total*((1/b)-(1/a)), strict=F)))*x
          # print(paste('f(',x,')=',num,'/',den,'=',num/den)) # = num*(1/den)*x
        }, 
        0, total
      )$value # else will return a list with value, abs.error, subdivisions, ...
    ) # i.e., integrate(PrSignalGivenTotal, /* from */ 0, /* to */ total);
    if(class(res) == 'try-error') {
      if(total > d*b) return(total-(d*b))
      else return(minx)
    } else {
      return(res)
    }
  }

} # }}}

## like gamma.mix, but stupider
gamma.ctl <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('methylated','unmethylated'), parallel=F){ # {{{

  if(parallel) require(multicore)
  if(is.null(channel)) { # {{{
    if(parallel) lstply <- mclapply else lstply <- lapply
    perchannel <- lstply(channels, function(channel) 
      gamma.ctl(object, channel=channel, allele=allele, parallel=parallel))
    names(perchannel) <- channels
    return(perchannel)
  } # }}}
  if(is.null(allele)) { # {{{
    if(parallel) lstply <- mclapply else lstply <- lapply
    perallele <- lstply(alleles, function(allele) {
      gamma.ctl(object, channel=channel, allele=allele, parallel=parallel)
    })
    names(perallele) <- alleles
    return(perallele) 
  } # }}}

  if(parallel) lstply <- mclapply else lstply <- lapply
  fg.params <- data.matrix(as.data.frame(lstply(1:dim(object)[2],function(i)
    gamma.mle(intensitiesByChannel(object[,i], channel, allele)))))
  bg.params <- gamma.bg(object, channel)
  colnames(fg.params) <- sampleNames(object)
  stopifnot(identical(colnames(fg.params), colnames(bg.params)))
  rownames(bg.params) <- c('bg.shape','bg.scale')
  rownames(fg.params) <- c('fg.shape','fg.scale')
  params <- t(rbind(fg.params, bg.params))

  signal <- data.matrix(as.data.frame(lstply(1:dim(object)[2], function(i) {
    sapply( intensitiesByChannel(object[,i], channel, allele), function(x) {
      gamma.conditional(x, params[i, ])
    })
  })))
  colnames(signal) <- sampleNames(object)
  rownames(signal) <- featureNames(object)[getProbesByChannel(object,channel)]
  return(signal)

} # }}}

## FIXME: add a qa step for the remapped beta-mixture scheme or don't use it
gamma.mix <- gamma.signal <- function(object, channel=NULL, channels=c('Cy3','Cy5'), parallel=F) { # {{{

## FIXME: use the lstply() hack everywhere else, too :-)
  if(parallel) require(multicore)
  if(is.null(channel)) { # {{{
    if(parallel) lstply <- mclapply else lstply <- lapply
    perchannel <- lstply(channels, function(channel) 
      gamma.mix(object, channel=channel, parallel=parallel))
    names(perchannel) <- channels
    return(perchannel)
  } # }}}

  both.params <- gamma.fg(object, channel)
  fg.params <- both.params$signal
  bg.params <- both.params$noise
  stopifnot(identical(colnames(fg.params), colnames(bg.params)))
  rownames(bg.params) <- c('bg.shape','bg.scale')
  rownames(fg.params) <- c('fg.shape','fg.scale')
  params <- t(rbind(fg.params, bg.params))
  ch <- channel
    
  ## here is where we diverge from gamma convolution against bg controls...
  ## instead of intensitiesByChannel, we have to use pseudo-totals from allelic
  ints <- intensitiesByChannel(object, channel)

  ## FIXME: use pvec or just straight C++
  if(parallel) lstply <- mclapply else lstply <- lapply
  signal <- lstply( names(ints), function(allele) {
    persubject <- data.matrix(as.data.frame(lstply(1:dim(object)[2],function(i){
      sapply( ints[[allele]][, i], function(x) {
        gamma.conditional(x, params[i, ]) ## vectorize!
      })
    })))
    colnames(persubject) <- sampleNames(object)
    rownames(persubject) <- featureNames(object)[getProbesByChannel(object,ch)]
    persubject
  })
  names(signal) <- names(ints)
  return(signal)

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

## FIXME: adjust negative and positive controls along with analytic probes
gamma.bgcorr <- function(object, how='mixture', offset=15, parallel=F) { # {{{
  
  if(annotation(object)=='HumanMethylation450k') stop('450ks not supported yet')
  else history.submitted <- as.character(Sys.time())

  ## FIXME: use switch()
  if( how %in% c('mixture','mix') ) {
    message('Control probes are not currently adjusted in the mixture model')
    signal <- gamma.mix(object, parallel=parallel)
    history.command <- "Applied gamma mixture model background correction."
  } else if( how %in% c('controls','ctl') ) {
    signal <- gamma.ctl(object, parallel=parallel)
    history.command <- "Applied gamma negative control background correction."
  }
  Ms <- methylated(object)
  Us <- unmethylated(object)
  for( ch in names(signal) ) {
    Ms[ getProbesByChannel(object, ch), ] <- signal[[ch]][['methylated']]
    Us[ getProbesByChannel(object, ch), ] <- signal[[ch]][['unmethylated']]
  }
  cloned <- clone(object)
  methylated(cloned) <- Ms + offset
  unmethylated(cloned) <- Us + offset
  betas(cloned) <- methylated(cloned) / total.intensity(cloned)
  pval.detect(cloned) <- 0.05 # resets betas, or at least, it should
  history.finished <- t.finish()
  cloned@history<- rbind(cloned@history,
                         data.frame(submitted=history.submitted,
                                    finished=history.finished,
                                    command=history.command))
  return(cloned)

} # }}}
