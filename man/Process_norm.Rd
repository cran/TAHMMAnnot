\name{Process_norm}
\alias{Process_norm}
\title{Normalization procedure}
\description{
  The whole normalization process from the pair files to the files with 3 columns, with ANOVA normalization or not, and with average on the dye-swap or not.
}
\usage{
  Process_norm(nbchip,nbarray,chip1_array1_Int1,chip1_array1_Int2,
chip1_array2_Int1,chip1_array2_Int2,chip1_array3_Int1,chip1_array3_Int2,
chip2_array1_Int1,chip2_array1_Int2,chip2_array2_Int1,chip2_array2_Int2,
chip2_array3_Int1,chip2_array3_Int2,fileOUT,fileOUTtext,fileOUTdata,
graph,Norm,MeanDyeswap,name.chr,sep.read="\t")

}
\arguments{
  \item{nbchip              }{The number of chip: 2 if dye-swap, 1 otherwise. }
  \item{nbarray             }{The number of array per chip.  }
  \item{chip1_array1_Int1   }{The name of the corresponding pair file.  }
  \item{chip1_array1_Int2   }{The name of the corresponding pair file.  }
  \item{chip1_array2_Int1   }{The name of the corresponding pair file.  }
  \item{chip1_array2_Int2   }{The name of the corresponding pair file.  }
  \item{chip1_array3_Int1   }{The name of the corresponding pair file.  }
  \item{chip1_array3_Int2   }{The name of the corresponding pair file.  }
  \item{chip2_array1_Int1   }{The name of the corresponding pair file.  }
  \item{chip2_array1_Int2   }{The name of the corresponding pair file.  }
  \item{chip2_array2_Int1   }{The name of the corresponding pair file.  }
  \item{chip2_array2_Int2   }{The name of the corresponding pair file.  }
  \item{chip2_array3_Int1   }{The name of the corresponding pair file.  }
  \item{chip2_array3_Int2   }{The name of the corresponding pair file.  }
  \item{fileOUT             }{The name of the file formatted for the ANOVA normalization.  }
  \item{fileOUTtext         }{The name of the text file with the results of the normalization. }
  \item{fileOUTdata         }{The name of the output file with normalized data.  }
  \item{graph               }{Logical. If we want to get all the graphs.  }
  \item{Norm                }{Logical. If we want to perform the ANOVA normalization. }
  \item{MeanDyeswap         }{Logical. If we want to average the intensities on the dye swap. }
  \item{name.chr            }{A vector of the name of the chromosomes to analyse.  }
  \item{sep.read            }{The field separator string to use in the \link{read.table} function. }
}
\details{
  Call the function Normalisation\_ANOVA if the option is TRUE
}
\value{
  The files formatted for the analyses with 3 columns, by chromosome and by strand.
}
\references{ 

 }
\author{ C. Berard, caroline.berard@agroparistech.fr    }
\seealso{Normalisation\_ANOVA}

\keyword{ }
\keyword{ }
