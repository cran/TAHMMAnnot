\name{ForwardAnnotR}
\alias{ForwardAnnotR}
\title{Forward step of the Forward-Backward algorithm in the case of model with annotation}
\description{
  Call the C function doing the Forward step of the Forward-Backward algorithm in the case of model with annotation  
}
\usage{
  ForwardAnnotR(Phi.vect,muHMM.vect,Mat.trans.norm.vect,annot,n,K,p)
}
\arguments{
  \item{Phi.vect       }{Vector. The bidimensionnal Gaussian emission densities. }
  \item{muHMM.vect    }{Vector. The stationary distribution of the heterogeneous Markov chain transformed in vector by column.  }
  \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{annot}{Vector of annotation.}
  \item{n}{The number of observations.}
  \item{K}{The number of groups.}
  \item{p}{Number of annotation categories.}
}
\details{
  
}
\value{
  F vector. Useful to calculate the posterior probabilities in the backward step.
}
\references{ 

 }

\seealso{BackwardAnnotR,ForwardR,ForwardAnnotrestart0R}

\keyword{internal}

