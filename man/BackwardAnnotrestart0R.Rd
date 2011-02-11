\name{BackwardAnnotrestart0R}
\alias{BackwardAnnotrestart0R}
\title{Backward step of the Forward-Backward algorithm in the case of model with annotation and with restart 0 option for the HMM.}
\description{
  Call the C function doing the Backward step of the Forward-Backward algorithm in the case of model with annotation and with restart 0 option for the HMM
}
\usage{
  BackwardAnnotrestart0R(muHMM.vect,F,Mat.trans.norm.vect,annot,diffannot,n,K,p)
}
\arguments{
  \item{muHMM.vect}{The stationary distribution of the heterogeneous Markov chain transformed in vector by column.}
  \item{F}{Vector. Output value of the ForwardR function.}
  \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{annot}{Vector of annotation. Each annotation category is transformed in an integer between 1 and p.}
  \item{diffannot}{Vector of iterated differences of annot. Useful to know the position where the annotation changes.}
  \item{n}{The number of observations.}
  \item{K}{The number of groups.}
  \item{p}{Number of annotation categories.}
}
\details{
  
}
\value{
  - G vector.
  - tau vector corresponds to the posterior probabilities vector.
}
\references{ 

 }

\seealso{ForwardAnnotrestart0R,BackwardAnnotR,BackwardR}

\keyword{internal}

