\name{gammaMME}
\alias{gammaMME}
\title{
  gammaMME: get the MME of a gamma distribution via fast conditional likelihood
}
\description{
  Method of moments estimator for the parameters of a gamma distribution
}
\usage{
  gammaMME(x)
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
  The method of moments estimate of the parameters.
}
\author{
  Tim Triche, Jr..
}
\examples{
  
  foo <- rgamma(100, 12, 10)
  gammaMME(foo)

}
