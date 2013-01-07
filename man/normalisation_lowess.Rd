\name{normalisation_lowess}
\alias{normalisation_lowess}
\title{Lowess Normalization}
\description{
  Application of a normalization by Lowess to remove technical biases. Formatting files from pair files to files with 3 columns, with Lowess normalization or not, and with average on the dye swap or not.
}
\usage{
  normalisation_lowess(fileIN,MeanDyeswap=TRUE,Lowess=TRUE,graph=TRUE,name.chr,underscore=FALSE,sep.write="\t",sep.read="\t")

}
\arguments{
 \item{fileIN        }{The file "lame.txt" which contains the number of the array (or the name of the pair file). The first number in the file is those corresponding to the reference sample in green (Cy3).  }
 \item{MeanDyeswap   }{Logical. If TRUE, the intensities are averaged on the dye swap.  }
 \item{Lowess        }{Logical. If TRUE, perform the Lowess normalization.  }
 \item{graph         }{Logical. If TRUE, print graphs before and after normalisation by chip, by array and by treatment for a dye (default value = FALSE).  }
 \item{name.chr      }{A vector of the name of the chromosomes to analyse.  }
 \item{underscore    }{Logical. If TRUE, there is an underscore in the probe ID.   } 
 \item{sep.write     }{The string to use for decimal points in the \link{write.table} function. }
 \item{sep.read      }{The field separator string to use in the \link{read.table} function.}

}
\details{
  WARNING: The models with annotation require marking information and work only if the repeat probes are removed.
}
\value{
  The files formatted for the analyses with 3 columns, by chromosome and by strand.
}


\seealso{Normalisation\_ANOVA, Process\_norm}
\author{ C. Berard, caroline.berard@agroparistech.fr    }

