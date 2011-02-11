\name{tau}
\alias{tau}
\title{Computation of the posterior probabilities}
\description{
  Computation of the posterior probabilities for a bidimensionnal Gaussian emission density
}
\usage{
  tau(x,sigma,mu,Pi)
}
\arguments{
  \item{x    }{Matrix (n*2) of data.  }
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

\seealso{tauannot}

\keyword{internal}
