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
    
double NormaliseF0(double *Mu,double **Phi, int nbClass)  
{
	//Rprintf("NormaliseF0\n");
    //Déclaration et initialisation des variables locales;   
    int y;
    double tot;
    tot=0;
    
    //Boucle permettant de calculer la somme des valeurs prises par F0 pour les K classes;
    for(y=0; y<nbClass; y++){ 
             tot=tot+Mu[y]*Phi[0][y];  
    }
    
   //Retourne la somme des éléments;  
   return tot;
}

/******************************************************************************/
/*                    Calcul du F0 normalisé                                  */
/******************************************************************************/        
    
void CalculF0(double *Mu,double **Phi,double **F, int nbClass)  
{
	//Rprintf("CalculF0\n");
   //Déclaration et initialisation des variables locales;   
   int y;
   double NormF0 = 0;
   
   //Appel de la fonction NormaliseF0 pour le calcul de la constante de normalisation;
   NormF0 =  NormaliseF0(Mu,Phi,nbClass);
   
   //Boucle permettant de calculer F0. normalisé;
   for(y=0; y<nbClass; y++){  
            F[0][y] =  Mu[y]*Phi[0][y]/NormF0;  
   }

}


/******************************************************************************/
/*            Calcul de la somme du produit entre Pi et F                     */
/******************************************************************************/


double sumFPi(double **Pi,double **F, int debut, int fin,int j,int t)  
{
	//Rprintf("sumFPi\n");
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;

   //Boucle permettant de calculer la somme produit entre Pi et F;
   for(x=debut; x<fin; x++){ 
             res=res+(Pi[x][j])*(F[t][x]);
   }

   //Retourne la somme des éléments;
   return res;
}    
   

/******************************************************************************/
/*             Calcul de la constante de normalisation pour F                 */
/******************************************************************************/    
    
double Normalise(double **Pi,double **Phi,double **F, int nbClass, int times)  
{
	//Rprintf("Normalise\n");
    //Déclaration et initialisation des variables locales;   
    int y;
    double tot;
    tot=0;
    
    //Boucle permettant de calculer la somme des valeurs prises par F pour les K classes;
    for(y=0; y<nbClass; y++){ 
          tot = tot + Phi[times][y]*sumFPi(Pi,F,0,nbClass,y,times-1);
    }
    
    //Retourne la somme des éléments;   
    return tot;
}

/******************************************************************************/
/*                     Calcul du F normalisé                                  */
/******************************************************************************/  
   
     void CalculF(double **Pi,double **Phi,double **F, int nbClass,int nbInd)  
{
	//Rprintf("CalculF\n");
     //Déclaration et initialisation des variables locales;   
     int x,y;
     double NormF = 0;
     
     //Double boucle permettant de calculer F à tout les instants t(sauf le premier) et pour toutes les classes;
     for(x=1; x<nbInd; x++){
              //Appel de la fonction Normalise pour le calcul de la constante de normalisation;
              NormF =  Normalise(Pi,Phi,F,nbClass,x); 
                     for(y=0; y<nbClass; y++){        
                             F[x][y] =  Phi[x][y]*sumFPi(Pi,F,0,nbClass,y,x-1)/NormF;    
                    }         
     }
}

/******************************************************************************/
/*       Etape Forward de l'algorithme Forward-Backward pour les HMM          */
/******************************************************************************/
    
        void Forward(double *PhiVect, double *Mu, double *PiVect,int *nbInd, int *nbClass, double *FVect)  
{
	//Rprintf("Forward\n");
	  int i;
          double **Phi = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                    Phi[i] = (double*)malloc(*nbClass*sizeof(double));
            }
    
    
          double **F = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   F[i] = (double*)malloc(*nbClass*sizeof(double));
            }
            
          double **Pi = (double**)malloc(*nbClass*sizeof(double*));
             for(i=0;i<*nbClass;i++){ 
                   Pi[i] = (double*)malloc(*nbClass*sizeof(double));
             }    
              
    VectToMat(PhiVect,Phi, *nbClass,*nbInd);
    VectToMat(PiVect,Pi, *nbClass,*nbClass);
    VectToMat(FVect,F, *nbClass,*nbInd);
   
    
    // Calcul du premier terme des F;
    CalculF0(Mu,Phi,F,*nbClass);
   
    // Calcul des autres termes(de l'instant 2 à t);           
    CalculF(Pi,Phi,F,*nbClass,*nbInd);
 
   
   MatToVect(Phi,PhiVect,*nbClass,*nbInd);
   MatToVect(F,FVect,*nbClass,*nbInd);
   MatToVect(Pi,PiVect,*nbClass,*nbClass);

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

    


}
  


