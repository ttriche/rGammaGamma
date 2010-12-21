plot.intensities <- function(object, samples=NULL, bybase=T) { # {{{
  colrs <- c(A='red',C='green',T='red')
  if(is.null(samples)) samples <- 1:dim(object)[2]
  bases <- sapply(names(colrs), function(b) which(fData(object)$Next_Base==b))
  meanb <- rowMeans(betas(object)[,samples], na.rm=T)
  pct <- list(low=unlist(which(meanb<0.5)), high=unlist(which(meanb>=0.5)))
  par(mfrow=c(ifelse(bybase,3,1),4))

  makeplots <- function(probes, pct, b='(any)', colr=NULL) {
    for(avg in names(pct)) {
      if( b=='(any)' ) colr <- ifelse(avg=='low','blue','orange')
      hist( methylated(object)[probes %i% pct[[avg]], samples],
            main=paste(avg, 's, base',b), xlim=c(0,2**16),
            breaks=50, col=colr, xlab='intensity')
    }
    for(avg in names(pct)) {
      if( b=='(any)' ) colr <- ifelse(avg=='low','blue','orange')
      hist(unmethylated(object)[probes %i% pct[[avg]], samples],
           main=paste(avg, 'beta Us, base',b), xlim=c(0,2**16),
           breaks=50, col=colr, xlab='intensity')
    }
  }

  if(bybase==T) 
    for(b in names(bases)) makeplots(bases[[b]], pct, b, colrs[[b]])
  else 
    makeplots(1:dim(object)[1], pct)
  
} # }}}

off.intensities <- function(object, samples=NULL, cutpoint=0.5, plots=F) { # {{{

  if(is.null(samples)) samples <- 1:dim(object)[2]
  colrs <- lapply(levels(fData(object)$Color_Channel), function(l) 
                  which(fData(object)$Color_Channel==l))
  names(colrs) <- levels(fData(object)$Color_Channel)
  U <- apply(betas(object)[,samples], 2, function(x) which(x < cutpoint))
  M <- apply(betas(object)[,samples], 2, function(x) which(x >= cutpoint))
  negs <- list(Red=Cy5.negs(object)[,samples], Grn=Cy3.negs(object)[,samples])
  norms <- list(Red=Cy5.norm(object)[,samples], Grn=Cy3.norm(object)[,samples])
  castoffs <- lapply(colrs, function(chan) {
    sapply(seq_along(samples), function(s) {
      c( methylated(object)[ U[[s]] %i% chan , samples[s] ],
         unmethylated(object)[ M[[s]] %i% chan , samples[s] ] )
    })
  })
  population <- lapply(colrs, function(chan) {
    sapply(seq_along(samples), function(s) {
      c( methylated(object)[ M[[s]] %i% chan , samples[s] ],
         unmethylated(object)[ U[[s]] %i% chan , samples[s] ] )
    })
  })

  # Now fit a gamma-normal convolution to each of the above (?)
  # It seems obvious to me that a mixture model with priors could do better
  subsets <- list(population=population,castoffs=castoffs,negs=negs,norms=norms)
  
  return(subsets) # red castoffs, green castoffs, red genpop, green genpop

} # }}}

# the so-called 'error function' (Abramowitz and Stegun 29.2.29)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

# and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)

# the canonical gamma MME -- turns out the MLE is actually fast and tractable!
gamma.mme <- function(x) { # {{{
  return(c(shape=(mean(as.matrix(x),na.rm=T)/sd(as.matrix(x),na.rm=T))^2,
           scale=var(as.matrix(x),na.rm=T)/mean(as.matrix(x),na.rm=T)))
} # }}}

# very fast approximation to the full Gamma MLE via Minka (2002) at MS Research
gamma.mle <- function(x,niter=100,tol=0.000000001) { # {{{

  meanlogx <- mean(log(na.omit(x)))
  meanx <- mean(x, na.rm=T)
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
  return(apply(negs(object, channel), 2, function(z) gamma.mle(z)))
} # }}}

gamma.fg <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5'), alleles=c('methylated','unmethylated')) { # {{{
  if(is.null(channel)) {
    perchannel <- lapply(channels, function(x) gamma.fg(object, x))
    names(perchannel) <- channels
    return(perchannel)
  }
  if(is.null(allele)) {
    perallele <- lapply(alleles, function(x) gamma.fg(object, channel, x))
    names(perallele) <- alleles
    return(perallele)
  }
  fn <- match.fun(paste(channel,allele,sep='.'))
  return(apply(fn(object), 2, function(z) gamma.mle(pmax(z,1))))
} # }}}

normal.bg <- function(object, channel=NULL, channels=c('Cy3','Cy5')) { # {{{
  if(is.null(channel)) { 
    perchannel <- lapply(channels, function(x) normal.bg(object,x))
    names(perchannel) <- channels
    return(perchannel)
  }
  return(colMeans(negs(object, channel), na.rm=T))
} # }}}

gamma.allelic <- function(object, channel=NULL, channels=c('Cy3','Cy5')) { # {{{
  if(is.null(channel)) {
    perchannel <- lapply(channels, function(x) gamma.fg(object, x))
    names(perchannel) <- channels
    return(perchannel)
  }
  stop('allelicity has not been correctly implemented yet')

  fn <- match.fun(paste(channel,allele,sep='.'))
  return(apply(fn(object), 2, function(z) gamma.mle(pmax(z,1))))
} # }}}

gamma.gamma.signal <- function(object, channel=NULL, allele=NULL, channels=c('Cy3','Cy5')) { # {{{
  stop('Gamma-Gamma convolution is not fully implemented yet')
  if(is.null(channel)) lapply(channels, function(x) gamma.bg.mode(object,x))
} # }}}

gamma.normal.params <- function(tot, bg, robust=F, mme=F) { # {{{ 

  if(is.null(tot)||is.null(bg)) stop('Non-NULL total and bg intensities needed')
  if(anyMissing(tot) || min(tot) < 1) {
    require(impute)
    tot <- pmax(impute.knn(tot)$data, 1)
  }
  if(class(tot) %in% c('matrix','data.frame')) {
    res <- t(sapply(colnames(tot), function(x) {
      gamma.normal.params(tot[,x],bg[,x],robust=robust)
    }))
    res <- cbind(res, colMeans(bg, na.rm=T), colSds(tot, na.rm=T))
    colnames(res) <- c('shape','scale','bg.mean','bg.sd')
    modes <- apply(res, 1, gamma.mode)
    res <- cbind(res, modes)
    colnames(res)[ncol(res)] <- 'mode'
    return(res)
  }
  ests <- c(mean(bg, na.rm=T), sd(bg, na.rm=T))
  if(robust) ests <- unlist(huber(na.omit(bg)))
  m <- m0 <- mu <- ests[1]
  s <- s0 <- ests[2]
  if( mme ) est <- gamma.mme(c(tot,bg))
  else est <- gamma.mle(c(tot,bg))
  a <- a0 <- est[1]
  b <- b0 <- est[2]

  # could maximize the joint loglikelihood for j negative controls and i probes
  params <- sapply(c(shape=a, scale=b, fg.mean=a*b, bg.mean=m, bg.sd=s), 
                   as.numeric)
  names(params) <- c('shape','scale','fg.mean','bg.mean','bg.sd')
  return(params)

} # }}}

gamma.gamma.params <- function(tot, bg, robust=F, mme=F) { # {{{ 

  if(is.null(tot)||is.null(bg)) stop('Non-NULL total and bg intensities needed')
  if(anyMissing(tot) || min(tot) < 1) {
    require(impute)
    tot <- pmax(impute.knn(tot)$data, 1)
  }
  if(class(tot) %in% c('matrix','data.frame')) {
    res <- t(sapply(colnames(tot), function(x) {
      gamma.normal.params(tot[,x],bg[,x],robust=robust)
    }))
    res <- cbind(res, colMeans(bg, na.rm=T), colSds(tot, na.rm=T))
    colnames(res) <- c('shape','scale','bg.mean','bg.sd')
    modes <- apply(res, 1, gamma.mode)
    res <- cbind(res, modes)
    colnames(res)[ncol(res)] <- 'mode'
    return(res)
  }
  ests <- c(mean(bg, na.rm=T), sd(bg, na.rm=T))
  if(robust) ests <- unlist(huber(na.omit(bg)))
  m <- m0 <- mu <- ests[1]
  s <- s0 <- ests[2]
  if( mme ) est <- gamma.mme(c(tot,bg))
  else est <- gamma.mle(c(tot,bg))
  a <- a0 <- est[1]
  b <- b0 <- est[2]

  # could maximize the joint loglikelihood for j negative controls and i probes
  params <- sapply(c(shape=a, scale=b, fg.mean=a*b, bg.mean=m, bg.sd=s), 
                   as.numeric)
  names(params) <- c('shape','scale','fg.mean','bg.mean','bg.sd')
  return(params)

} # }}}

allelicity <- function(x, betacut=0.5) { # {{{ mcut=0, derp
  allele <- apply(betas(x), 2, function(y) y >= betacut)
  return(list(mostlysignal=(allele*methylated(x))+((1-allele)*unmethylated(x)),
              mostlynoise=((1-allele)*methylated(x))+(allele*unmethylated(x))))
} # }}}

# do this once per channel per allelicity (RedMostlySignal, RedMostlyNoise, ...)
gamma.normal.signal <- function (x, params)  { # {{{

  if(any(params<=0)) stop('All parameters must be positive')
  bg.mean <- params['bg.mean']
  bg.sd <- params['bg.sd']
  shape <- params['shape']
  scale <- params['scale']

  # could use the mode below in a gamma-gamma model
  mu.sf <- x - bg.mean - (bg.sd**2)*scale
  expectedSgivenT <- function(x, params) {
    # ... this is where the approximation comes in ... 
  }

  signal <- expectedSgivenT(x, params)

  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal

} # }}}

# normal-exponential deconvolution (WEHI code) for comparison with the above
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

# this is the final result of all the above chicanery... a normal-gamma model
mcsv.normgamma <- function(object,channel=NULL,allele=NULL,channels=c('Cy3','Cy5'),alleles=c('mostlysignal','mostlynoise'),annotation=NULL,betacut=0.5,robust=F, use.score=T, parallel=F, offset=50){ #{{{

  x <- clone(object)
  if(any(methylated(x)<=0)) methylated(x)[which(methylated(x)==0)] <- 1
  if(any(unmethylated(x)<=0)) unmethylated(x)[which(unmethylated(x)==0)] <- 1
  betas(x) <- methylated(x) / (methylated(x) + unmethylated(x))
  if( is.null(QCdata(x)) ) # can adapt this to crlmm no problem...
    stop('Your object does not have any control bead data! Exiting.')

  history.submitted <- as.character(Sys.time())

  # this can and should be a methylumi or Lumi method
  if(is.null(annotation)) annotation <- annotation(x)
  if(is.null(annotation)) annotation <- 'IlluminaHumanMethylation27k'
  require(paste(annotation,'db',sep='.'), character.only=T) # for COLORCHANNEL
  annochan <- paste(annotation,'COLORCHANNEL',sep='')
  color <- unlist(mget(featureNames(x), get(annochan), ifnotfound=NA))
  cy5 <- which(color=='Red')
  cy3 <- which(color=='Grn')
  probes <- list(Cy5=cy5, Cy3=cy3, cy5=cy5, cy3=cy3) # backwards compatibility
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

# for comparison:

