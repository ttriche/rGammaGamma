require('gsl') # for hyperg_1F1 in the conditional expectation of the signal
require('methylumIDAT') # for the various utility functions, e.g. cy3(), cy5()

################################################################################

# the canonical gamma MME -- but the cMLE is lower variance and just as fast!
gamma.mme <- function(x) { # {{{
  return(c(shape=(mean(as.matrix(x),na.rm=T)/sd(as.matrix(x),na.rm=T))^2,
           scale=var(as.matrix(x),na.rm=T)/mean(as.matrix(x),na.rm=T)))
} # }}}

# very fast approximation to the full Gamma MLE via Minka (2002) at MS Research
gamma.mle <- function(x,niter=100,tol=0.000000001,minx=1) { # {{{

  meanlogx <- mean(log(na.omit(pmax(x,minx))))
  meanx <- mean(pmax(x,minx), na.rm=T)
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

## FIXME:
##
## Get rid of the dependence on methylumiCSV's Cy3.negs and Cy5.negs methods
##
gamma.bg <- function(object, channel=NULL, channels=c('Cy3','Cy5')) { # {{{
  if(is.null(channel)) {
    perchannel <- lapply(channels, function(x) gamma.bg(object,x))
    names(perchannel) <- channels
    return(perchannel)
  }
  return(apply(negs(object, channel), 2, function(z) gamma.mle(z)))
} # }}}

gamma.bg.ebayes <- function(object, channel=NULL, channels=c('Cy3','Cy5')){ #{{{
  
  # pool the gamma parameter distributions to estimate hyperparameters 
  
  # integrate(marginalScale, 0, Inf)
  # integrate(marginalShape, 0, Inf)

} # }}}

allelic <- allelicity <- function(x, channel=NULL, allele=NULL, hard=F) { # {{{ 

  # we could fit a fuzzy mixture model here, w/binary pi0
  if(!is.null(channel) && tolower(channel) %in% c('cy3','cy5')) {
    getchan <- match.fun(tolower(channel))
    probes <- getchan(x)
  } else {
    probes <- 1:dim(x)[1]
  }

  # if we fit the fuzzy mixture, add a pi1 attribute
  pi0 <- methylated(x)[probes,] / total.intensity(x)[probes,]
  # also, if we do that, should we weight by detection here?

  if(hard) pi0 <- round(pi0) # hard == dichotomized 

  intensities <-list(
    signal=((pi0*methylated(x)[probes,]) + ((1-pi0)*unmethylated(x)[probes,])),
     noise=(((1-pi0)*methylated(x)[probes,]) + (pi0*unmethylated(x)[probes,]))
  )
  
  if(is.null(allele)) return(intensities) 
  else return(intensities[[allele]])

} # }}}

gamma.allelic <- gamma.fg <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('signal','noise'), hard=F) { # {{{

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

  apply(allelic(object,channel=channel,allele=allele,hard=hard),2,gamma.mle)

} # }}}

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

gamma.signal <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('signal','noise'), hard=F) { # {{{

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
  sapply(1:dim(totals)[2], function(j) {
    sapply(totals[-pars,j], gamma.conditional, params=totals[pars,j])
  })

  ## Now reconstruct the methylated and unmethylated bg-corrected signals,
  ## using the adjusted betas to push them towards the "appropriate" place.

} # }}}

# gamma-gamma deconvolution (my code) for comparison; mvalcut == log2it(betacut)
mcsv.gammagamma <- function(object,channel=NULL,allele=NULL,channels=c('Cy3','Cy5'),alleles=c('mostlysignal','mostlynoise'),annotation=NULL,betacut=0.5,robust=F, use.score=T, parallel=F, offset=25){ #{{{

  x <- clone(object)
  methylated(x) <- pmax(methylated(x), 1)
  unmethylated(x) <- pmax(unmethylated(x), 1)
  betas(x) <- methylated(x) / (methylated(x) + unmethylated(x))
  if(is(x, 'MethyLumiM')) qcdata <- controlData(x)
  if(is(x, 'MethyLumiSet')) qcdata <- QCdata(x)
  if( is.null(QCdata(x)) ) {
    stop('Your object does not appear to have any control bead data! Exiting.')
  } else {
    history.submitted <- as.character(Sys.time())
  }

  # this can and should be a methylumi or Lumi method
  if(is.null(annotation)) annotation <- annotation(x)
  if(is.null(annotation)) annotation <- 'IlluminaHumanMethylation27k'
  require(paste(annotation,'db',sep='.'), character.only=T)
  probes <- list(Cy5=cy5(x), Cy3=cy3(x))

  if(ncol(betas(x)) != ncol(negs(x, 'Cy3'))) {
    stop("Control beads do not match the number of samples in the dataset!!!")
  } else {
    cat("Background mean and sd estimated from", nrow(negs(x,'Cy3')),"probes\n")
  }

  # should use correlation structure here, or at least merge the probes;
  # when we regress m.u.corr ~ gc.content (weighted by bisulfite conversion), 
  # it explains a massive amount of the differences observed...
  #
  # for allelic deconvolution, this means ( M<0.5, U>0.5 ) and ( M>0.5, U<0.5 )
  # works fine with beta.logit(betas(object))>0 or exprs(object)>logit(betacut)
  #
  allelic <- allelicity(x, betacut)
  allelic.chan <- lapply(channels, function(ch) {
    by.allele <- lapply(alleles, function(allele) {
      list(tot=allelic[[allele]][probes[[ch]],],bg=negs(x,ch),fg=nonnegs(x,ch))
    })
    names(by.allele) <- alleles
    return(by.allele)
  })
  names(allelic.chan) <- channels
  if(parallel) {
    require(multicore)
    options("cores"=2) # trust me here...
    params <- mclapply(allelic.chan, function(ch) {
      mclapply(ch, function(al.ch) {
        gamma.normal.params(
          al.ch[['tot']],al.ch[['bg']],robust=robust,use.score=use.score
        )
      })
    })
  } else { 
    params <- lapply(allelic.chan, function(ch) {
      lapply(ch, function(al.ch) {
        gamma.normal.params(
          al.ch[['tot']],al.ch[['bg']],robust=robust,use.score=use.score
        )
      })
    })

  }
  return(params)

  stop('normal-gamma deconvolution is not yet finished due to an 1F1 problem!')

  #############################################################################
  ## THIS IS AS FAR AS I HAVE GOTTEN IN TERMS OF ADAPTING NORMEXP TO NORMGAM ##
  #############################################################################

  m.u.chan <- F # mapping: merged probes to unmerged probes, sample X channel 

  zz <- zmu <- zsigma <- zalpha <- list()
  for( channel in c('cy3','cy5') ) {
    dat <- rbind(dat.cpg[[channel]][['m']], dat.cpg[[channel]][['u']])
    m <- nrow(dat) 
    n <- ncol(dat)
    mu <- sigma <- alpha <- rep(NA, n)
    z <- rbind(dat, dat.pctrl[[channel]], dat.nctrl[[channel]])
    for(i in 1:n){
      #if(robust){ 

      #} else { 
        mu[i] <- mean(dat.nctrl[[channel]][, i])
        sigma[i] <- sd(dat.nctrl[[channel]][, i])
        alpha[i] <- max(mean(dat[, i]) - mu[i], 10)
      #}
      z[, i] <- mcsv.nexp.signal(c(mu[i], log(sigma[i]), log(alpha[i])), z[,i])
    }
    z <- z + offset
    zz[[channel]] <- z
    zmu[[channel]] <- mu
    zsigma[[channel]] <- sigma
    zalpha[[channel]] <- alpha
  }

  # jointly estimated
  m.x <- methylated(x)
  u.x <- unmethylated(x)
  for( channel in c('cy3','cy5') ) {
    m.x[ probes[[channel]], ] <- zz[[channel]][mprobes[[channel]], ]
    u.x[ probes[[channel]], ] <- zz[[channel]][uprobes[[channel]], ]
    dat.ctrl[[channel]][pos, ] <- zz[[channel]][pcprobes[[channel]], ]
    dat.ctrl[[channel]][neg, ] <- zz[[channel]][ncprobes[[channel]], ]
  }
  methylated(x) <- m.x
  unmethylated(x) <- u.x
  dat.nctrl <- list(cy5=Cy5(ctrl)[neg,], cy3=Cy3(ctrl)[neg,])
  dat.pctrl <- list(cy5=Cy5(ctrl)[pos,], cy3=Cy3(ctrl)[pos,])
  pData(x) <- cbind(pData(x), mu.Cy5=zmu$cy5, mu.Cy3=zmu$cy3, 
                              sigma.Cy5=zsigma$cy5, sigma.Cy3=zsigma$cy3,
                              alpha.Cy5=zalpha$cy5, alpha.Cy3=zalpha$cy3)
  Cy3(ctrl) <- dat.ctrl$cy3
  Cy5(ctrl) <- dat.ctrl$cy5
  QCdata(x) <- ctrl # adjusted 

	varMetadata(x)[rownames(varMetadata(x)) == "mu.Cy5", ] <- 
    "Estimated mean background intensity, red channel"
	varMetadata(x)[rownames(varMetadata(x)) == "mu.Cy3", ] <- 
    "Estimated mean background intensity, green channel"
	varMetadata(x)[rownames(varMetadata(x)) == "sigma.Cy5", ] <- 
    "Estimated standard deviation of background intensities, red channel"
	varMetadata(x)[rownames(varMetadata(x)) == "sigma.Cy3", ] <- 
    "Estimated standard deviation of background intensities, green channel"
	varMetadata(x)[rownames(varMetadata(x)) == "alpha.Cy5", ] <- 
    "Estimated mean signal intensity, red channel"
	varMetadata(x)[rownames(varMetadata(x)) == "alpha.Cy3", ] <- 
    "Estimated mean signal intensity, green channel"

  history.command <- "Applied normexp background correction."
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)
  zz <- zmu <- zsigma <- zalpha <- list()
  for( channel in c('cy3','cy5') ) {
    dat <- rbind(dat.cpg[[channel]][['m']], dat.cpg[[channel]][['u']])
    m <- nrow(dat) 
    n <- ncol(dat)
    mu <- sigma <- alpha <- rep(NA, n)
    z <- rbind(dat, dat.pctrl[[channel]], dat.nctrl[[channel]])
    for(i in 1:n){
      #if(robust){ 

      #} else { 
        mu[i] <- mean(dat.nctrl[[channel]][, i])
        sigma[i] <- sd(dat.nctrl[[channel]][, i])
        alpha[i] <- max(mean(dat[, i]) - mu[i], 10)
      #}
      z[, i] <- mcsv.nexp.signal(c(mu[i], log(sigma[i]), log(alpha[i])), z[,i])
    }
    z <- z + offset
    zz[[channel]] <- z
    zmu[[channel]] <- mu
    zsigma[[channel]] <- sigma
    zalpha[[channel]] <- alpha
  }

} # }}}

################################################################################

# normal-exponential deconvolution (WEHI code) for comparison with the above
normal.bg <- function(object, channel=NULL, channels=c('Cy3','Cy5')) { # {{{
  if(is.null(channel)) { 
    perchannel <- lapply(channels, function(x) normal.bg(object,x))
    names(perchannel) <- channels
    return(perchannel)
  }
  return(colMeans(negs(object, channel), na.rm=T))
} # }}}

mcsv.nexp.signal <- function (par, x)  { # {{{
  mu <- par[1]
  sigma <- exp(par[2])
  sigma2 <- sigma * sigma
  alpha <- exp(par[3])
  if (alpha <= 0) 
    stop("alpha must be positive")
  if (sigma <= 0) 
    stop("sigma must be positive")
  mu.sf <- x - mu - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(dnorm(0, mean = mu.sf, sd = sigma, log = TRUE)-
    pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log = TRUE))
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal
} # }}}

mcsv.gamma.signal <- function (par, x)  { # {{{
  mu <- par[1]
  sigma <- exp(par[2])
  sigma2 <- sigma * sigma
  alpha <- exp(par[3])
  if (alpha <= 0) 
    stop("alpha must be positive")
  if (sigma <= 0) 
    stop("sigma must be positive")
  mu.sf <- x - mu - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(dnorm(0, mean = mu.sf, sd = sigma, log = TRUE)-
    pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log = TRUE))
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal
} # }}}
