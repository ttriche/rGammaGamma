\name{rGammaGamma-package}
\alias{rGammaGamma-package}
\alias{rGammaGamma}
\docType{package}
\title{
  Implementation of a gamma deconvolution for background correction.
}
\description{
  SNP and methylation microarrays exhibit a different distribution of signal
  and background intensities than expression microarrays.  The flexibility of
  the gamma distribution suits the observed signal and background intensity 
  distributions well, encompassing as a special case the normal-exponential 
  convolution (the de facto standard for expression background correction).
  In addition, a signal/noise mixture model is implemented for empirical 
  deconvolution of cross-hybridization between probe alleles.
}
\details{
\tabular{ll}{
Package:  rGammaGamma
Type:     Package
Version:  1.0.11
Date:     2012-12-04
License:  Artistic
LazyLoad: yes
}
}
\author{
Tim Triche, Jr.<ttriche@usc.edu>
}
\keyword{ package }
\examples{
## not run

##  foo <- methylumIDAT(barcodes, parallel=T)
##  foo.bg <- methylumi.bgcorr(foo, method='goob')
}
