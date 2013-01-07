\name{HMM}
\alias{HMM}
\title{ Hidden Markov Model for tiling array analysis }
\description{
  A 4-state hidden Markov model with bidimensional Gaussian emission densities used to model the joint distribution of two signals where the dependence between neighboring probes is taken into account. Analysis for one chromosome. 
}
\usage{
HMM(color=c("navajowhite","grey","black","green","red"),
fileIN="MoyDye_CHR01.txt",var1="IS1",var2="IS2",K=4,eps=10e-6,
file.RData="Param_CHR01.RData",file.graph="Graph_CHR01.png",
file.flagdb="Fichier_flagdb.txt",file.OUT="Output_CHR01.txt",
file.INIT="Graph_init_CHR01.png",file.param="Fichier_Param_CHR01.txt",
threshold=FALSE,s=0.7,random.init=FALSE,a.inf=60,a.sup=75,b.inf=5,
b.sup=20,PtInit=8,theta=0.4,int.max=8,max.iter=1000,header=TRUE,
sep = "\t", ...)
}
\arguments{
  \item{color        }{A vector of K+1 colors to color the graph with the different groups.  }
  \item{fileIN       }{The name of the Input file with 3 columns: ID, var1, var2.  }
  \item{var1         }{The name of the first column of intensities.  }
  \item{var2         }{The name of the second column of intensities.  }
  \item{K            }{The number of components.  }
  \item{eps          }{The stopping relative precision limit for stopping EM algorithm.  }
  \item{file.RData   }{The name of the RData object to save the parameter estimates.  }
  \item{file.graph   }{The name of the graph.  }
  \item{file.flagdb  }{The name of the flagdb file in order to visualize the results in Flagdb++.  }
  \item{file.OUT     }{The name of the output file.  }
  \item{file.INIT    }{The name of the graph file with the initialized parameters.  }
  \item{file.param   }{The name of the text file with estimator values.  }
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
  The Input file "FileIN" has to be ordered according to the position of the probes along the chromosome. 
  The parameters are estimated using the EM algorithm. The E step is achieved with the Forward-Backward algorithm which is implemented in C.
  
}
\value{
   A RData object with the parameter estimations. 
   A graph png with the visualization of the initialization. 
   A graph png to visualize the probes after classification. 
   An output file with the posterior probabilities and the declared status for each probe using the classification rule.
   A text file with the parameter estimations.
   A flagdb file which has the specific form to be input in the Flagdb++ software.
}
\references{ 
    Unsupervised Classification for Tiling Arrays: ChIP-chip and Transcriptome. Caroline Berard, Marie-Laure Martin-Magniette, Veronique Brunaud, Sebastien Aubourg, Stephane Robin, Statistical Applications in Genetics and Molecular Biology, 10(50), 2011.
 }


\author{ C. Berard, caroline.berard@agroparistech.fr }
\seealso{HMM\_Annot, Mixture\_Annot, Bidim\_Gaussian\_Mixture, HMM\_Annot\_restart0}

