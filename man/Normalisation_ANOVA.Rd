\name{Normalisation_ANOVA}
\alias{Normalisation_ANOVA}
\title{ANOVA Normalization}
\description{
  Application of a normalization by ANOVA to remove technical biases
}
\usage{
   Normalisation_ANOVA(fileIN,fileOUTtext,fileOUTdata,graph)
}
\arguments{
   \item{fileIN}{The name of the input data file with 6 columns named ID, CHIP, ARRAY, DYE, TREATMENT, VALUE. }
   \item{fileOUTtext}{The name of the results file for estimation parameters and means squares. }
   \item{fileOUTdata}{The name of the output data file with 6 columns (the first 5 columns are the same as in input file and the column 6 is named VALUEnorm which is the intensity of the probe after normalisation). }
   \item{graph}{Logical. If TRUE, print graphs before and after normalisation by chip, by array and by treatment for a dye (default value = FALSE). }
}
\details{
  It must be at least 2 chips or 2 arrays and 2 dyes.
  This function is called in the Process\_norm.R function.
}
\value{
  - One file for the parameters estimation of the normalization function.
  - One file for the output data after normalization.
  - If the option graph=TRUE, pdf graphs before and after normalization.
}

\author{ C. Berard, caroline.berard@agroparistech.fr    }

