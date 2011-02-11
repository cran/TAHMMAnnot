\name{HMM_allchr}
\alias{HMM_allchr}
\title{Hidden Markov Model for a set of chromosomes defined by the user.}
\description{
  Do the HMM function for a set of chromosomes defined by the user.
}
\usage{
   HMM_allchr(name.chr=c("CHR01","CHR02","CHR03","CHR04","CHR05"),
color=c("navajowhite","grey","black","green","red"),
beginfileIN="MoyDye",var1="IS1",var2="IS2",K=4,eps=10e-6,
beginfile.RData="Param",beginfile.graph="Graph",
file.flagdb="Fichier_flagdb.txt",beginfile.OUT="Output",
beginfile.INIT="Graph_init",beginfile.param="Fichier_Param",
threshold=FALSE,s=0.7,random.init=FALSE,a.inf=60,a.sup=75,
b.inf=5,b.sup=20,PtInit=8,theta=0.4,int.max=8,max.iter=1000,
header=TRUE,sep = "\t", ...)

}
\arguments{
  \item{name.chr          }{A vector of the name of the chromosomes to analyse.  }
  \item{color             }{A vector of K+1 colors to color the graph with the different groups.  }
  \item{beginfileIN       }{The beginning of the name of the Input file with 3 columns which is added the name of the chromosome.  }
  \item{var1              }{The name of the first column of intensities.  }
  \item{var2              }{The name of the second column of intensities.  }
  \item{K                 }{The number of components.  }
  \item{eps               }{The stopping relative precision limit for stopping EM algorithm.  }
  \item{beginfile.RData   }{The beginning of the name of the RData object to save the parameter estimates which is added the name of the chromosome.  }
  \item{beginfile.graph   }{The beginning of the name of the graph file which is added the name of the chromosome.  }
  \item{file.flagdb       }{The name of the flagdb file in order to visualize the results in Flagdb++.  }
  \item{beginfile.OUT     }{The beginning of the name of the output file which is added the name of the chromosome.  }
  \item{beginfile.INIT    }{The beginning of the name of the graph file with the initialized parameters which is added the name of the chromosome.  }
  \item{beginfile.param   }{The beginning of the name of the text file with values of parameters which is added the name of the chromosome.  }
  \item{threshold    }{Logical. If we need a classification threshold to classify with the MAP rule. }
  \item{s            }{The value of the classification threshold.  }
  \item{random.init  }{Logical. If we want a random initialization or not.  }
  \item{a.inf        }{The minimum value of the proportion of the noise group in order to initialize.  }
  \item{a.sup        }{The maximum value of the proportion of the noise group in order to initialize.  }
  \item{b.inf        }{The minimum value of the proportion of the identical group in order to initialize.  }
  \item{b.sup        }{The maximum value of the proportion of the identical group in order to initialize.  }
  \item{PtInit       }{The coordinate of intersection of lines in order to initialize. }
  \item{theta        }{The value of the angle between lines in order to initialize. The higher the angle is, more the lines are close.  }
  \item{int.max      }{The maximum intensity value of the noise group.  }
  \item{max.iter     }{The maximum number of iterations to do before stopping the EM algorithm.  }
  \item{header       }{ a logical value indicating whether the file contains the names of the variables as its first line. Used in the \link{read.table} function. }
  \item{sep          }{ the field separator string to use in the \link{read.table} function. }
  \item{\dots        }{ Further arguments to be passed to the \link{read.table} function. }

}
\details{
  This function just call the HMM function for each chromosome defined in the name.chr variable.
}
\value{
  The values are those of the HMM function for each chromosome.
}
\references{ 
Unsupervised Classification for Tiling Arrays: ChIP-chip and Transcriptome. Caroline Berard, Marie-Laure Martin-Magniette, Veronique Brunaud, Sebastien Aubourg, Stephane Robin, arXiv:1104.5429v1.

 }

\author{ C. Berard, caroline.berard@agroparistech.fr }
\seealso{HMM\_Annot\_allchr, Melange\_Annot\_allchr, Melange\_gaussien\_bidim\_allchr, \cr
HMM\_Annot\_restart0\_allchr}

\keyword{ }
\keyword{ }
