\name{ForwardR}
\alias{ForwardR}
\title{Forward step of the Forward-Backward algorithm }
\description{
  Call the C function doing the Forward step of the Forward-Backward algorithm 
}
\usage{
  ForwardR(Phi.vect,muHMM,Mat.trans.norm.vect,n,K)
}
\arguments{
  \item{Phi.vect       }{Vector. The bidimensionnal Gaussian emission densities. }
  \item{muHMM       }{Vector. The stationary distribution of the Markov chain.  }
  \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{n}{The number of observations.}
  \item{K}{The number of groups.}
}

\value{
  F vector. Useful to calculate the posterior probabilities in the backward step.
}

\seealso{BackwardR,ForwardAnnotR,ForwardAnnotrestart0R}

\keyword{internal}

