\name{gamma.mle}
\alias{gamma.mle}
\title{
  gamma.mle: get the MLE of a gamma distribution via fast conditional likelihood
}
\usage{
  gamma.mle(x)
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
  gamma.mle(foo)

}
