\name{gammaIntegral}
\alias{gammaIntegral}
\title{
  gammaIntegral: get the conditional expectation of a convolved gamma RV
}
\description{
  Conditional expectation of one (out of two) gamma RVs given their sum & params
}
\usage{
  gammaIntegral(total, params, offset=50, minx=1)
}
\arguments{
  \item{total}{
    a vector of sums (convolutions)
  }
  \item{params}{
    parameters for the signal and background distributions
  }
  \item{offset}{
    an offset to be added to the estimate
  }
  \item{minx}{
    a minimum value below which we will replace x 
  }
}
\details{
  See Triche et al., Nucleic Acids Research (submitted).
}
\value{
  The conditional expectation(s) of the signal given the total.
}
\author{
  Tim Triche, Jr..
}
