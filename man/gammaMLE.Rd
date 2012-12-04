\name{gammaMLE}
\alias{gammaMLE}
\title{
  gammaMLE: get the MLE of a gamma distribution via fast conditional likelihood
}
\description{
  Maximum likelihood estimator for the parameters of a gamma distribution
}
\usage{
  gammaMLE(x, w=NULL, niter=100, tol=0.1, minx=1)
}
\arguments{
  \item{x}{
    A vector 
  }
  \item{w}{
    Weights
  }
  \item{niter}{
    Maximum number of iterations
  }
  \item{tol}{
    Maximum difference in parameters at convergence
  }
  \item{minx}{
    Minimum permissible value for x (smaller values will be bumped to this)
  }
}
\details{
  Not yet.
}
\value{
  The MLE of the parameters.
}
\author{
  Tim Triche, Jr..
}
\examples{
  
  foo <- rgamma(100, 12, 10)
  gammaMLE(foo, w=NULL, niter=100, tol=0.1, minx=1)

}
