\name{tauannot}
\alias{tauannot}
\title{Computation of the posterior probabilities in the case of a model with annotation}
\description{
  Computation of the posterior probabilities for a bidimensionnal Gaussian emission density in the case of a model with annotation
}
\usage{
  tauannot(xbis,sigma,mu,Pi)
}
\arguments{
  \item{x    }{Matrix (n*3). The two first colums are data. The third column is the annotation.  }
  \item{sigma}{Variance matrix of the Gaussian.}
  \item{mu   }{Mean vector of the Gaussian.  }
  \item{Pi   }{Vector of the proportion in each group. }
}
\details{
  
}
\value{
  A n*K matrix of posterior probabilities.
}
\references{ 

 }

\seealso{tau}

\keyword{internal}
