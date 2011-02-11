\name{BackwardAnnotR}
\alias{BackwardAnnotR}
\title{Backward step of the Forward-Backward algorithm in the case of model with annotation}
\description{
  Call the C function doing the Backward step of the Forward-Backward algorithm in the case of model with annotation  
}
\usage{
  BackwardAnnotR(F,Mat.trans.norm.vect,annot,n,K,p)
}
\arguments{
 
  \item{F}{Vector. Output value of the ForwardR function.}
  \item{Mat.trans.norm.vect   }{Vector. The transition matrix of the Markov chain transformed in vector by column.  }
  \item{annot}{Vector of annotation.}
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

\seealso{ForwardAnnotR,BackwardR,BackwardAnnotrestart0R}

\keyword{internal}

