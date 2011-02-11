\name{BackwardR}
\alias{BackwardR}
\title{Backward step of the Forward-Backward algorithm }
\description{
  Call the C function doing the Backward step of the Forward-Backward algorithm 
}
\usage{
  BackwardR(F,Mat.trans.norm.vect,n,K)
}
\arguments{
  \item{F}{Vector. Output value of the ForwardR function.}
  \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{n}{The number of observations.}
  \item{K}{The number of groups.}
}
\details{
  
}
\value{
  - G vector.
  - tau vector corresponds to the posterior probabilities vector.
}
\references{ 

 }

\seealso{ForwardR,BackwardAnnotR,BackwardAnnotrestart0R}

\keyword{internal}

