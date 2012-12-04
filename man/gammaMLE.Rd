\name{gammaMLE}
\alias{gammaMLE}
\title{
  gammaMLE: get the MLE of a gamma distribution via fast conditional likelihood
}
\description{
  Maximum likelihood estimator for the parameters of a gamma distribution
}
\usage{
  gammaMLE(x)
}
\arguments{
  \item{x}{
    A vector 
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
  gammaMLE(foo)

}
