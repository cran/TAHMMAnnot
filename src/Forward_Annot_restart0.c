#include <stdio.h>  /* directives au préprocesseur */
//#include <iostream>
#include <math.h>
#include <stdlib.h>
//#include <fstream>
//#include <string>
#include <Rinternals.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

#include <R.h>
#include "Convert.h"



/******************************************************************************/
/*             Calcul de la constante de normalisation pour F0                */
/******************************************************************************/    
    
double NormaliseF0Annot_restart0(double **Mu,double **Phi, int *annot, int nbClass)  
{
	//Rprintf("NormaliseF0\n");
    //Déclaration et initialisation des variables locales;   
    int y;
    double tot;
    tot=0;
    
    //Boucle permettant de calculer la somme des valeurs prises par F0 pour les K classes;
    for(y=0; y<nbClass; y++){ 
	 
             tot=tot+Mu[y][annot[0]-1]*Phi[0][y];  
    }
    
   //Retourne la somme des éléments;  
   return tot;
}

/******************************************************************************/
/*                    Calcul du F0 normalisé                                  */
/******************************************************************************/        
    
void CalculF0Annot_restart0(double **Mu,double **Phi,double **F, int *annot, int nbClass)  
{
	//Rprintf("CalculF0\n");
   //Déclaration et initialisation des variables locales;   
   int y;
   double NormF0 = 0;
   
   //Appel de la fonction NormaliseF0 pour le calcul de la constante de normalisation;
   NormF0 =  NormaliseF0Annot_restart0(Mu,Phi,annot,nbClass);
   
   //Boucle permettant de calculer F0. normalisé;
   for(y=0; y<nbClass; y++){  
            F[0][y] =  Mu[y][annot[0]-1]*Phi[0][y]/NormF0;  
   }

}


/******************************************************************************/
/*            Calcul de la somme du produit entre Pi et F                     */
/******************************************************************************/


double sumFPiAnnot_restart0(double **Mu,double **Pi,double **F, int *annot, int *diffannot, int debut, int fin,int j,int t)  
{
	//Rprintf("sumFPi\n");
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;
   //int cpt=0;
   //Boucle permettant de calculer la somme produit entre Pi et F;
   for(x=debut; x<fin; x++){ 
         //     cpt =cpt+1;
         //     if (t<6)
         //     {
	//	printf("%d\t", cpt);
	//	printf("%d\t", t);
	//	printf("%d\t", annot[t+1]);
	//	printf("%d\t", diffannot[t]);
	//      printf("%f\t", Pi[x][(annot[t+1]-1)*fin+j]); 
	//	printf("%f\t", F[t][x]);
	//	printf("%f\t", (Pi[x][(annot[t+1]-1)*fin+j]*(F[t][x]))); 
	//	printf("%f\n", Mu[x][annot[t+1]-1]);
        //      }
		if(diffannot[t] == 0)
		{
                   res=res+(Pi[x][(annot[t+1]-1)*fin+j])*(F[t][x]);
 		}
		else
		{
 		   res=res+(Mu[x][annot[t+1]-1])*(F[t][x]);
		}
   }

   //Retourne la somme des éléments;
   return res;
}    
  

/******************************************************************************/
/*             Calcul de la constante de normalisation pour F                 */
/******************************************************************************/    
    
double NormaliseAnnot_restart0(double **Mu,double **Pi,double **Phi,double **F, int *annot,int *diffannot,int nbClass, int times)  
{
	//Rprintf("Normalise\n");
    //Déclaration et initialisation des variables locales;   
    int y;
    double tot;
    tot=0;
    
    //Boucle permettant de calculer la somme des valeurs prises par F pour les K classes;
    for(y=0; y<nbClass; y++){ 
          tot = tot + Phi[times][y]*sumFPiAnnot_restart0(Mu,Pi,F,annot,diffannot,0,nbClass,y,times-1);
    }
    
    //Retourne la somme des éléments;   
    return tot;
}

/******************************************************************************/
/*                     Calcul du F normalisé                                  */
/******************************************************************************/  
   
     void CalculFAnnot_restart0(double **Mu,double **Pi,double **Phi,double **F, int *annot, int *diffannot, int nbClass,int nbInd)  
{
	//Rprintf("CalculF\n");
     //Déclaration et initialisation des variables locales;   
     int x,y;
     double NormF = 0;
     
     //Double boucle permettant de calculer F à tout les instants t(sauf le premier) et pour toutes les classes;
     for(x=1; x<nbInd; x++){
              //Appel de la fonction Normalise pour le calcul de la constante de normalisation;
              NormF =  NormaliseAnnot_restart0(Mu,Pi,Phi,F,annot,diffannot,nbClass,x); 
                     for(y=0; y<nbClass; y++){        
                             F[x][y] =  Phi[x][y]*sumFPiAnnot_restart0(Mu,Pi,F,annot,diffannot,0,nbClass,y,x-1)/NormF;    
                    }         
     }
}

/******************************************************************************/
/*       Etape Forward de l'algorithme Forward-Backward pour les HMM          */
/******************************************************************************/
    
        void Forward_Annot_restart0(double *PhiVect, double *MuVect, double *PiVect,int *annot,int *diffannot,int *nbInd, int *nbClass, int *p, double *FVect)  
{
	//Rprintf("Forward\n");
	  int i;
          double **Phi = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                    Phi[i] = (double*)malloc(*nbClass*sizeof(double));
            }
    
	  double **Mu = (double**)malloc(*nbClass*sizeof(double*));
            for(i=0;i<*nbClass;i++){ 
                    Mu[i] = (double*)malloc(*p*sizeof(double));
            }
    
          double **F = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   F[i] = (double*)malloc(*nbClass*sizeof(double));
            }
            
          double **Pi = (double**)malloc(*nbClass*sizeof(double*));
             for(i=0;i<*nbClass;i++){ 
                   Pi[i] = (double*)malloc((*nbClass*(*p))*sizeof(double));
             }    
              
    VectToMat(PhiVect,Phi, *nbClass,*nbInd);
    VectToMat(PiVect,Pi, *nbClass*(*p),*nbClass);
    VectToMat(FVect,F, *nbClass,*nbInd);
    VectToMat(MuVect,Mu, *p, *nbClass);

    //rprintf("These are two other numbers: %d %d", k, 4);
 //int j;
   // for (i=0;i<*nbClass;i++)
   // {
   // for (j=0;j<*nbClass*3;j++)
   //  printf("%f\t", Pi[i][j]);
   //  printf("\n" ); 
  //  }

//for (i=0;i<*nbClass;i++)
  //  {
    // for (j=0;j<3;j++)
   //  printf("%f\t", Mu[i][j]);
   //  printf("\n" ); 
  //  }

    // Calcul du premier terme des F;
      CalculF0Annot_restart0(Mu,Phi,F,annot,*nbClass);
   
    // Calcul des autres termes(de l'instant 2 à t);           
      CalculFAnnot_restart0(Mu,Pi,Phi,F,annot,diffannot,*nbClass,*nbInd);
 
   
   MatToVect(Phi,PhiVect,*nbClass,*nbInd);
   MatToVect(F,FVect,*nbClass,*nbInd);
   MatToVect(Pi,PiVect,*nbClass*(*p),*nbClass);
   MatToVect(Mu,MuVect,*p,*nbClass);

   //libérer l'espace mémoire
     for(i=0;i<*nbInd;i++){ 
         free(Phi[i]) ;
     }
    free(Phi);
    
    for(i=0;i<*nbClass;i++){ 
         free(Pi[i]) ;
    }
    free(Pi);
 
    for(i=0;i<*nbInd;i++){ 
         free(F[i]) ;
    }
    free(F);

    for(i=0;i<*nbClass;i++){ 
         free(Mu[i]) ;
    }
    free(Mu);

}
  


