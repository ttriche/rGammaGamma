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

## faster approximation (in C++, but it seems to screw up multicore somehow)
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

## for plotting and comparing the results of gamma.mixest vs. real data
make.fn <- function(ests,channel,allele,bgorfg,type=c('rgamma','dgamma')){ # {{{
  if(type=='rgamma') {
    return(
      function(n) {
        rgamma(n, shape=ests[paste(channel,allele,bgorfg,'shape',sep='.')],
                  scale=ests[paste(channel,allele,bgorfg,'scale',sep='.')])
      }
    )
  }
  if(type=='dgamma') {
    return(
      function(x) {
        rgamma(x, shape=ests[paste(channel,allele,bgorfg,'shape',sep='.')],
                  scale=ests[paste(channel,allele,bgorfg,'scale',sep='.')])
      }
    )
  }
} # }}}
make.dgamma.fn <- function(ests,channel,allele,bgorfg) { # {{{
  make.fn(ests, channel, allele, bgorfg, type='dgamma')
} # }}}
make.rgamma.fn <- function(ests,channel,allele,bgorfg) { # {{{
  make.fn(ests, channel, allele, bgorfg, type='rgamma')
} # }}}

## FIXME: definitely move this to C++!
gamma.mixest <- function(object, ch=NULL, al=NULL, sp=T, w.beta=F, als=c('methylated','unmethylated'), chs=c('Cy3','Cy5'), cuts=c(.15,.85),use.CpGi=T, confine=F, verbose=F ){ # {{{
  if(is.null(ch)) {
    perchannel <- lapply(chs, function(x) gamma.nonspecific(object,x,al,sp=sp))
    return(cbind(perchannel[[1]], perchannel[[2]]))
  }
  if(sp && is.null(al)) {
    perallele <- lapply(als, function(a) gamma.nonspecific(object,ch,a,sp=sp))
    return(cbind(perallele[[1]], perallele[[2]])) 
  }
  if(ch=='Cy3') probes <- cy3(object)
  if(ch=='Cy5') probes <- cy5(object)
  ann <- paste(annotation(object),'ISCPGISLAND',sep='')
  if(use.CpGi) {
    CpGi <- which(unlist(mget(featureNames(object),get(ann),ifnotfound=NA))==1)
    nonCpGi<-which(unlist(mget(featureNames(object),get(ann),ifnotfound=NA))==0)
  } else { 
    CpGi <- nonCpGi <- 1:dim(object)[1]
  }
  other <- list(methylated='unmethylated', unmethylated='methylated')
  alnm <- lapply(names(other), function(al) toupper(substr(al, 1, 1)))
  names(alnm) <- names(other)

  # this next part might be ideal for multicore
  nonspecific <- sapply(1:dim(object)[2], function(i) {
    probes <- list(unmethylated=(which(betas(object)[,i] < cuts[1]) %i% CpGi),
                   methylated=(which(betas(object)[,i] > cuts[2]) %i% nonCpGi))
    if( w.beta ) {
      piM <- (1-betas(object)[low,i])
      piU <- betas(object)[high,i]
    } else {
      piM <- piU <- 1
    }
    pi0 <- list(methylated=piM, unmethylated=piU)
    if(sp) { # split M and U 

      ## FIXME: should I use > 0.85 and < 0.15 to estimate FG, or BG?!?
      if(confine) {
        ndist <- assayDataElement(object,other[[al]])[probes[[al]],i]
        sdist <- assayDataElement(object,al)[probes[[al]],i]
      } else { 
        ndist <- assayDataElement(object,other[[al]])[probes[[other[[al]]]],i]
        sdist <- assayDataElement(object,al)[probes[[al]],i]
      }
      if(verbose) {
        message('used ',length(ndist),' probes to estimate ',ch,' ',other[[al]],
                ' nonspecificity on ', sampleNames(object)[i])
        message('used ',length(sdist),' probes to estimate ',ch,' ', al, 
                ' specificity on ', sampleNames(object)[i])
      }
      c( gamma.cmle(ndist), gamma.cmle(sdist) )
    } else {
      ndist <- c(methylated(object)[low,i]*piM,unmethylated(object)[high,i]*piU)
      gamma.cmle(sdist)
    }
  })
  if(sp) {
    rownms <- c( paste(ch, alnm[[other[[al]]]],'bg',c('shape','scale'),sep='.'),
                 paste(ch, alnm[[al]], 'fg', c('shape','scale'), sep='.') )
    rownames(nonspecific) <- rownms # includes specific/fg params, too
  } else { 
    rownames(nonspecific) <- paste(ch, 'bg', c('shape','scale'), sep='.')
  }
  nonspec <- (t(nonspecific))
  rownames(nonspec) <- sampleNames(object)
  return(nonspec)
} # }}}

## FIXME: probably farm this out, too 
gamma.from.SE <- function(mu, s) { # {{{ bead-level model
  return( c(shape=(mu/s)*(mu/s), scale=(s*s)/mu) )
} # }}}

## FIXME: hook into methylumIDAT::pval.detect()
gamma.detect <- function(object, use.SE=F) { # {{{

  ## if estimates for the nonspecific intensity aren't already there, add them
  if(!all(c('Cy5.bg.shape','Cy5.bg.scale','Cy3.bg.shape','Cy3.bg.scale') %in%
          varLabels(object))) {
    pData(object) <- cbind(pData(object), gamma.nonspecific(object))
  }
  chs <- c('Cy3','Cy5')
  chfuns <- c('Cy3','Cy5')

  ## if we have standard errors and bead numbers, use the evidence ratio 
  if(all(c('methylated.SE','unmethylated.SE','methylated.N','unmethylated.N') 
         %in% assayDataElementNames(object))) { 
    message('Bead numbers and standard errors found, computing Pr(Xs>Xb)...') 
    stop('LRT not finished yet')
    pvals.by.channel <- lapply(chs, function(ch) {

    })
  } else { ## otherwise just compute the exceedence probability  
    message('Standard errors and bead numbers missing, using pgamma(mean)...') 
    pvals.by.channel <- lapply(chs, function(ch) {

      apply(pmax(methylated(object)[pch[[ch]],],
                 unmethylated(object)[pch[[ch]],]), c(1,2), function(x) {
            pgamma(x, shape=pshape[[ch]], scale=pscale[[ch]])
      })
    })
  }

  pvals.merged <- betas(object) # get the right size of matrix
  stop('pvalues now must be mapped back via cy3() and cy5()!')
  return(pvals.merged)

} # }}}

## FIXME: move this to C++
gamma.miller.pars <- function(vec) { # {{{
  c( p=prod(vec), q=sum(vec), r=length(vec), s=length(vec) )
} # }}}

## FIXME: move this to C++
gamma.miller.posterior <- function(g.and.a, params) { # {{{

  g0 = g.and.a[1]
  a0 = g.and.a[2]
  p = params[1]
  q = params[2]
  r = params[3]
  s = params[4]
  Z = gamma.miller.normconst(g0, p, q, r, s) 
  g = ((p**(g0-1))*gamma((s*g0)+1)) / (Z*(q**(1-(s*g0)))*(gamma(g0)**r))
  scalefn = NULL ## FINISH THIS
  a = integrate(scalefn, 0, Inf)
  return( shape=g, scale=a )

} # }}}

## FIXME: move this to C++
gamma.miller.normconst <- function(g, p, q, r, s) { # {{{
  fxn = function(g) ((p**(g-1))*gamma((s*g)+1))/((gamma(g)**r)*(q**((s*g)+1)))
  integrate(fxn, 0, Inf)  
} # }}}

## FIXME: perhaps do this in C++ for pvalues
gamma.x.gt.y <- function(x, params.x, params.y) { # {{{ thanks to John Cook
  g <- params.x[1] # fg shape
  a <- params.x[2] # fg scale
  d <- params.y[1] # bg shape
  b <- params.y[2] # bg scale
  pbeta( a/(a+b), d, g)
} # }}}

## FIXME: handle CpG stratification using (updated) annotations!
gamma.bg.ebayes <- function(object, channel=NULL, channels=c('Cy3','Cy5')){ #{{{

  if(is.null(channel)) { # {{{
    perchannel = lapply(channels, function(ch) gamma.bg.ebayes(object, ch))
    names(perchannel) <- channels
    return(perchannel)
  } # }}}

  # estimate priors from negative controls as outlined in Appendix I 
  priors = apply(negctls(object, channel), 2, gamma.miller.pars)
  params = apply(negctls(object, channel), 2, gamma.mle)
  both = cbind(priors, params)

  # normConstants = apply(both, 2, function(x) integrate(fZ, 0, Inf, params=x))
  # integrate(marginalScale, 0, Inf) # update the shape within each stratum
  # integrate(marginalShape, 0, Inf) # update the scale within each stratum

  ## return updated per-stratum estimates (I figure this goes in the QC pData) 
  return( posterior.params.per.property.plenum )

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

## just get the damned off-channel estimates
gamma.nonspecific <- function(object,ch=NULL,chs=c('Cy3','Cy5'),use.U=F,cuts=c(0.05, 0.85),CpG=F) { # {{{
  if(is.null(ch)) { # {{{
    perchannel <- lapply(chs, function(ch) {
      nonspec = gamma.nonspecific(object,ch=ch,use.U=use.U,cuts=cuts)
      colnames(nonspec) = paste(ch, colnames(nonspec), sep='.')
      nonspec
    })
    return(cbind(perchannel[[1]], perchannel[[2]]))
  } # }}}
  if(is.null(ch)) { # {{{
    perchannel <- lapply(chs, function(ch) {
      nonspec = gamma.nonspecific(object,ch=ch,use.U=use.U,cuts=cuts)
      colnames(nonspec) = paste(ch, colnames(nonspec), sep='.')
      nonspec
    })
    return(cbind(perchannel[[1]], perchannel[[2]]))
  } # }}}
  if(ch=='Cy3') probes <- cy3(object)
  if(ch=='Cy5') probes <- cy5(object)
  als = c( 'M', 'U' )
  pars = c( 'shape', 'scale' )
  nonspecific <- sapply(1:dim(object)[2], function(i) {
    bv = (pmax(methylated(object),1)/pmax(total.intensity(object),2))[,i]
    low = intersect(which( bv < cuts[1] ), probes)
    high = intersect(which( bv > cuts[2] ), probes)
    est = c(rep(gamma.mle(methylated(object[low,i])),2))
    if(use.U) est[3:4] = gamma.mle(unmethylated(object[high,i]))
    nms = unlist(sapply(als,function(x)paste(x,'bg',pars,sep='.'),simplify=F))
    names(est) = nms
    est
  })
  nonspec <- t(nonspecific)
  rownames(nonspec) = sampleNames(object)
  return(nonspec)
} # }}}

## simple: in-channel estimates for a given allele
gamma.specific <- function(object, ch=NULL, chs=c('Cy3','Cy5'), al=NULL, alleles=c('methylated','unmethylated')) { # {{{

  if(is.null(ch)) { # {{{
    perchannel <- lapply(chs, function(ch) {
      spec = gamma.specific(object, ch=ch, al=al)
      colnames(spec) = paste(ch, colnames(spec), sep='.')
      spec
    })
    return(cbind(perchannel[[1]], perchannel[[2]]))
  } # }}}
  if(is.null(al)) { # {{{
    perallele <- lapply(alleles, function(al) {
      spec = gamma.specific(object, ch=ch, al=al)
      colnames(spec) = paste(toupper(substr(al,1,1)), colnames(spec), sep='.')
      spec
    })
    return(cbind(perallele[[1]], perallele[[2]]))
  } # }}}
  if(ch=='Cy3') probes <- cy3(object)
  if(ch=='Cy5') probes <- cy5(object)
  pars = c( 'shape', 'scale' )
  specific <- sapply(1:dim(object)[2], function(i) {
    est = gamma.mle(assayDataElement(object, al)[probes,i])
    names(est) = paste('fg',pars,sep='.')
    est
  })
  spec <- t(specific)
  rownames(spec) = sampleNames(object)
  return(spec)
} # }}}

# combines the above two functions for easy cbind()'ing into pData
gamma.mixparams <- function(x, use.U=F) { # {{{
  cbind(gamma.nonspecific(x, use.U=use.U), gamma.specific(x))
} # }}}

## FIXME: move this to C++ as soon as humanly possible (ideally with matrix arg)
gamma.integral <- function(total, params, minx=1) { # {{{

  ## this bit is the most obvious "farm me out to C++" piece of all...
  if(length(total) > 1) return(sapply(total, gamma.integral, params=params))
  
  g = params[1]
  a = params[2]
  d = params[3]
  b = params[4]
  bg.mean = d * b 
  bg.sd = sqrt( d * b * b )
  cat('total =',total,'... bg.mean =',bg.mean,'... bg.sd =',bg.sd,"\n")
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
      return(pmax(total-bg.mean, minx))
    } else {
      return(pmax(res, minx))
    }
  }

} # }}}

## FIXME: add a qa step for the remapped beta-mixture scheme or don't use it
gamma.mix <- gamma.mix2 <- function(object, use.U=F, minx=15, parallel=F){ # {{{

  if( 'Cy3.fg.shape' %in% varLabels(object) ||  # {{{ don't do this twice...
      'M.Cy3.fg.shape' %in% varLabels(object) ) { 
    message('You seem to already have had background correction. Exiting.')
    return(object)
  } else { 
    pData(object) = cbind(pData(object), gamma.mixparams(object, use.U=use.U))
    params = list( 
      Cy3 = list(
        methylated = pData(object)[, c('Cy3.M.fg.shape','Cy3.M.fg.scale',
                                       'Cy3.M.bg.shape','Cy3.M.bg.scale')  ],
        unmethylated = pData(object)[, c('Cy3.U.fg.shape','Cy3.U.fg.scale',
                                         'Cy3.U.bg.shape','Cy3.U.bg.scale')]
      ),
      Cy5 = list(
        methylated = pData(object)[, c('Cy5.M.fg.shape','Cy5.M.fg.scale',
                                       'Cy5.M.bg.shape','Cy5.M.bg.scale')  ],
        unmethylated = pData(object)[, c('Cy5.U.fg.shape','Cy5.U.fg.scale',
                                         'Cy5.U.bg.shape','Cy5.U.bg.scale')]
      )
    ) 
    params = lapply(params, function(h) lapply(h, function(al) data.matrix(al)))
  } # }}}

  if(parallel) require(multicore)
  if(parallel) lstply = mclapply else lstply = lapply
  probes = list(Cy5=cy5(object), Cy3=cy3(object))
  alleles = list(M='methylated', U='unmethylated')
  channels = list(Cy5='Cy5', Cy3='Cy3') ## FIXME: add for 450k: , Both='New')
  signals = lstply( alleles, function( al ) {
    data.matrix(as.data.frame(lstply( 1:dim(object)[2], function( i ) { 
      scratch = assayDataElement(object, al)[ , i ]
      for( ch in channels ){
        chpr = probes[[ch]]
        scratch[chpr] = gamma.integral(scratch[chpr],params[[ch]][[al]][i,])
      }
      return(scratch)
    })))
  })

  x = clone(object)
  methylated(x) = signals$M
  unmethylated(x) = signals$U
  if(is(x, 'MethyLumiM')) {
    exprs(x) = log2(pmax(signals$M,1)/pmax(signals$U,1))
  } else if(is(x, 'MethyLumiSet')) {
    betas(x) = pmax(signals$M,1)/pmax(total.intensity(x),1)
  }
  return(x)

} # }}}

## FIXME: add a qa step for the remapped beta-mixture scheme or don't use it
gamma.mix1 <- function(object, channel=NULL, channels=c('Cy3','Cy5'), parallel=F, al=NULL, alleles=c('methylated','unmethylated'), sp=T, how='nonspecific') { # {{{

  ## FIXME: use the lstply() hack everywhere else, too :-)
  if(parallel) require(multicore)
  if(is.null(channel)) { # {{{
    if(parallel) lstply <- mclapply else lstply <- lapply
    perchannel <- lstply(channels, function(channel) 
      gamma.mix(object, channel=channel, parallel=parallel))
    names(perchannel) <- channels
    return(perchannel)
  } # }}}

  ## here is where we diverge from gamma convolution against bg controls...
  ## instead of intensitiesByChannel, we have to use pseudo-totals from allelic
  ints <- intensitiesByChannel(object, channel)

  ## FIXME: default to using gamma.nonspecific()
  if(sp) both.params <- gamma.nonspecific(object, channel)
  both.params <- gamma.fg(object, channel)
  fg.params <- both.params$signal
  bg.params <- gamma.nonspecific(object, channel)
  stopifnot(identical(colnames(fg.params), colnames(bg.params)))
  rownames(bg.params) <- c('bg.shape','bg.scale')
  rownames(fg.params) <- c('fg.shape','fg.scale')
  params <- t(rbind(fg.params, bg.params))
  ch <- channel

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

## FIXME: switch to using C++ and/or OpenMP to speed this up tolerably
## FIXME: adjust negative and positive controls along with analytic probes
## FIXME: add a log entry for gamma deconvolution and note how it was done
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
