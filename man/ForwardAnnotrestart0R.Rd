\name{ForwardAnnotrestart0R}
\alias{ForwardAnnotrestart0R}
\title{Forward step of the Forward-Backward algorithm in the case of model with annotation and with restart 0 option for the HMM.}
\description{
  Call the C function doing the Forward step of the Forward-Backward algorithm in the case of model with annotation and with restart 0 option for the HMM
}
\usage{
  ForwardAnnotrestart0R(Phi.vect,muHMM.vect,Mat.trans.norm.vect,annot,diffannot,n,K,p)
}
\arguments{
  \item{Phi.vect       }{Vector. The bidimensionnal Gaussian emission densities. }
  \item{muHMM.vect}{The stationary distribution of the heterogeneous Markov chain transformed in vector by column.}
   \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{annot}{Vector of annotation. Each annotation category is transformed in an integer between 1 and p.}
  \item{diffannot}{Vector of iterated differences of annot. Useful to know the position where the annotation changes.}
  \item{n}{The number of observations.}
  \item{K}{The number of groups.}
  \item{p}{Number of annotation categories.}
}

\value{
  F vector. Useful to calculate the posterior probabilities in the backward step.
}


\seealso{BackwardAnnotrestart0R,ForwardAnnotR,ForwardR}

\keyword{internal}

