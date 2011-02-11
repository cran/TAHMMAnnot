\name{matW}
\alias{matW}
\title{Estimator of the dispersion matrix of a group k. }
\description{
  Compute the product matrices for the estimator of the dispersion matrix of a group k.
}
\usage{
  matW(x,k,mu)
}
\arguments{
  \item{x    }{The matrix (n*2) of the data.  }
  \item{k    }{The number of the component.  }
  \item{mu   }{The mean vector of the bidimensionnal Gaussian density.  }
}
\details{
  We compute only the product of matrices. It is not exactly the estimator of W because the product with tau is missing.
}
\value{
  The product of matrices for the estimator of the dispersion matrix of a group k. 
}
\references{ 

 }

\seealso{}

\keyword{internal}
