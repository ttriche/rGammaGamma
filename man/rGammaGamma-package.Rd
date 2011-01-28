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
Version:  0.9
Date:     2011-01-27
License:  Artistic
LazyLoad: yes
}
}
\author{
Tim Triche, Jr.<ttriche@usc.edu>
}
\keyword{ package }
\examples{
  foo <- methylumIDAT(barcodes, parallel=T)
  foo.bg <- gammagamma(foo)
}
