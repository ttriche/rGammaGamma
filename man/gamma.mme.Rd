\name{gamma.mme}
\alias{gamma.mme}
\title{
  gamma.mme: get the MME of a gamma distribution via fast conditional likelihood
}
\usage{
  gamma.mme(x)
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
  gamma.mme(foo)

}
