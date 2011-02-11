#include <stdio.h>  /* directives au préprocesseur */
//#include <iostream>
#include <math.h>
#include <stdlib.h>
//#include <fstream>
//#include <string>
//#include <Rinternals.h>
//#include <gsl/gsl_sf_log.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_blas.h>

#include <R.h>
#include "Convert.h"


/******************************************************************************/
/*            Calcul de la somme intervenant dans la formule du Tau           */
/******************************************************************************/  

double CalculSommeTauAnnot_restart0(double **Mu,double **Pi,double **Tau,double **G, int *annot, int *diffannot, int debut, int fin,int j,int t)  
{
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;
   //int cpt=0;
     //Boucle permettant de calculer la somme produit entre Pi, Tau et 1/G;
     for(x=debut; x<fin; x++){ 
	//cpt =cpt+1;
        //      if (t>107026 & t<107032)
        //      {
	//	printf("%d\t", cpt);
	//	printf("%d\t", t);
	//	printf("%d\t", annot[t]);
	//	printf("%d\t", diffannot[t]);
	//      printf("%f\t", Pi[j][(annot[t]-1)*fin+x]); 
	//	printf("%f\t", Tau[t][x]); 
	//	printf("%f\t", (Pi[j][(annot[t]-1)*fin+x]*Tau[t][x])); 
	//	printf("%f\t", G[t][x]);
	//	printf("%f\t", (Pi[j][(annot[t]-1)*fin+x]*Tau[t][x]/G[t][x]));
	//	printf("%f\t", Mu[j][annot[t]-1]);
	//	printf("%f\t", (Mu[j][annot[t]-1]*Tau[t][x]));
	//	printf("%f\n", (Mu[j][annot[t]-1]*Tau[t][x]/G[t][x]));
          //    }

         if(diffannot[t-1] == 0)  //car length(diffannot) = length(annot)-1;
	 {
             res=res+((Pi[j][(annot[t]-1)*fin+x]*Tau[t][x])/G[t][x]);
 	 }
	 else
	 {
 	     res=res+((Mu[j][annot[t]-1]*Tau[t][x])/G[t][x]);
	 }
   }
   //Retourne la somme des éléments;
   return res;
}    


/******************************************************************************/
/*             Calcul de la somme du produit entre Pi et F,                   */
/*            identique à sumFPi sauf pour les indices du Pi                  */
/******************************************************************************/
    
  
  double sumFPi2Annot_restart0(double **Mu,double **Pi,double **F, int *annot, int *diffannot, int debut, int fin,int j,int t)  
{
   //Déclaration et initialisation des variables locales;      
   int x;
   double res;
   res=0;
   //int cpt=0;
   //Boucle permettant de calculer la somme produit entre Pi et F;
   for(x=debut; x<fin; x++){ 
   //cpt =cpt+1;
     //     if (t<6)
     //     {
	//	printf("%d\t", cpt);
	//	printf("%d\t", t);
	//	printf("%d\t", annot[t+1]);
	//	printf("%d\t", diffannot[t]);
	  //    printf("%f\t", Pi[x][(annot[t+1]-1)*fin+j]); 
	//      printf("%f\t", F[t][x]); 
	//	printf("%f\t", (Pi[x][(annot[t+1]-1)*fin+j])*(F[t][x])); 
	//	printf("%f\t", Mu[x][annot[t+1]-1]);
//		printf("%f\n", (Mu[x][annot[t+1]-1])*(F[t][x]));
 //          }
		
         if(diffannot[t] == 0)
	 {
             res=res+((Pi[x][(annot[t+1]-1)*fin+j])*(F[t][x]));
 	 }
	 else
	 {
 	     res=res+((Mu[x][annot[t+1]-1])*(F[t][x]));
	 }
   }
   //Retourne la somme des éléments;
   return res;
}    
  
/******************************************************************************/
/*                               Calcul du G                                  */
/******************************************************************************/  

void CalculGAnnot_restart0(double **Mu,double **Pi,double **F, int *annot, int *diffannot, int nbClass,int nbInd,double **G)  
{
     //Déclaration et initialisation des variables locales;   
     int x,y;
     
     //Double boucle permettant de calculer G à tout les instants et pour toutes les classes;
     for(x=0; x<(nbInd-1); x++){
                     for(y=0; y<nbClass; y++){                                      
                             G[x+1][y] =sumFPi2Annot_restart0(Mu,Pi,F,annot,diffannot,0,nbClass,y,x);   
			     G[0][y] = 0;                                
                    }         
     }   
}


/******************************************************************************/
/*                              Calcul des Tau                                */
/******************************************************************************/  

void CalculTauAnnot_restart0(double **Mu,double **Pi,double **F, int *annot, int *diffannot, int nbClass,int nbInd,double **G,double **Tau)  
{
     //Déclaration et initialisation des variables locales;   
     int x,y;
     double NormF = 0;
     
     //Calcul de la matrice G nécessaire au calcul des Tau;
     CalculGAnnot_restart0(Mu,Pi,F,annot,diffannot,nbClass,nbInd,G);
     
     //Double boucle permettant de calculer Tau à tout les instants t(sauf le dernier) et pour toutes les classes;     
     for(x=(nbInd-2); x>=0; x--){
                     for(y=0; y<nbClass; y++){                                   
                             Tau[x][y] =  F[x][y]*CalculSommeTauAnnot_restart0(Mu,Pi,Tau,G,annot,diffannot,0,nbClass,y,x+1);
                      }         
     }
}


/******************************************************************************/
/*         Calcul du premier instant de l'étape Backward  (Tau[N][.])         */
/******************************************************************************/

void CalculTauNAnnot_restart0(double **F,int nbInd, int nbClass, double **Tau)  
{
    int y;
    for(y=0; y<nbClass; y++){ 
    Tau[nbInd-1][y]=F[nbInd-1][y];
    }
}   


/******************************************************************************/
/*       Etape Backward de l'algorithme Forward-Backward pour les HMM         */
/******************************************************************************/
    
void Backward_Annot_restart0(double *MuVect,double *FVect, double *PiVect,int *annot,int *diffannot, int *nbInd, int *nbClass, int *p, double *TauVect, double *GVect)  
{
     int i;
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
             
        double **Tau = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   Tau[i] = (double*)malloc(*nbClass*sizeof(double));
  	    }
  
  double **G = (double**)malloc(*nbInd*sizeof(double*));
            for(i=0;i<*nbInd;i++){ 
                   G[i] = (double*)malloc(*nbClass*sizeof(double));
            }
  
      VectToMat(PiVect,Pi, *nbClass*(*p),*nbClass);
      VectToMat(FVect,F, *nbClass,*nbInd);    
      VectToMat(TauVect,Tau, *nbClass,*nbInd); 
      VectToMat(GVect,G, *nbClass,*nbInd); 
      VectToMat(MuVect,Mu, *p, *nbClass);

    //Premier instant de l'étape Backward (instant t); 
    CalculTauNAnnot_restart0(F,*nbInd,*nbClass,Tau); 
       
    //Calcul des autres termes(de l'instant 1 à t-1);           
    CalculTauAnnot_restart0(Mu,Pi,F,annot,diffannot,*nbClass,*nbInd,G,Tau);
   
    MatToVect(F,FVect,*nbClass,*nbInd);
    MatToVect(Pi,PiVect,*nbClass*(*p),*nbClass);
    MatToVect(Tau,TauVect,*nbClass,*nbInd);
    MatToVect(G,GVect,*nbClass,*nbInd);
    MatToVect(Mu,MuVect,*p,*nbClass);

 for(i=0;i<*nbInd;i++){ 
                   free(F[i]) ;
            }
    free(F);
 for(i=0;i<*nbInd;i++){ 
                   free(Tau[i]) ;
            }
    free(Tau);
 for(i=0;i<*nbInd;i++){ 
                   free(G[i]) ;
            }
    free(G);
 for(i=0;i<*nbClass;i++){ 
                   free(Pi[i]) ;
            }
    free(Pi);
  for(i=0;i<*nbClass;i++){ 
         free(Mu[i]) ;
    }
    free(Mu);


}


   
