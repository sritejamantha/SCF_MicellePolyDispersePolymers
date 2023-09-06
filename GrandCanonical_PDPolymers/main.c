/* SCF in 1D Spherical coordinates ---> Sriteja Mantha, Shuqnhu Qi, and Friedrike Schmid */
/* System ---> Polydisperse polymers with Nn=10.0 in a spherical solvent box of radius Lx */
/* Grand canonical ensemble */
/* Exponential of Sol Chem Pot set to 1.0*/
/* Exponential of Chem Pot corresponding to Pol with chain length N taken from corresponding, */
/* single chain SCF calculations in mixed ensemble at same PhiPbar as current GC simulations*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include<gsl/gsl_sf_gamma.h>

#include "systemF.h"
#include "getDensity.h"
#include "energy.h"
#include "dataOutput.h"

double PI=3.141592653589793;
clock_t start, finish;

int main()
{
 int i;
 double runtime;
 systemF s;
 int it;
 long int iter=100000000; /*iteration steps to solve SCF equations*/
 double errorA, errorB, errorS;
 double freeEnergy=0.0, freeEnergyNew=0.0,freeEnergyDiff=0.0;
 double accuracy=1.e-9;

 start=clock();

 s.ns=1000;             /*Chain grid points to solve diffusion equations*/
 s.f=0.50;              /*Fraction of solvophobic groups in the chain*/
 s.Nn=1;                /*Scaled Nn to compute Schulz-Zimm distribution*/

 s.Nscale=10.0;         /*Chain length scaling factor*/ 

 s.Lx=10.0;             /*radius of the sphere within which SCF is solved*/
 s.nx=100;              /*Grid points along the spatial direction*/

 s.lambda=0.01;         /*Lambda factor to obtain new potential values*/
  
 s.xAB=10.0;            /*Flory parameters, scaled to Nscale*/
 s.xAS=20.0;
 s.xBS=-0.50;
 s.kappa=100.0;         /*Compressibility constant, also scaled to Nscale*/

 s.vol=(4.0*PI*pow(s.Lx,3.0))/3.0;                /*sphere volume*/
 s.dV=s.vol/((4.0*PI*pow(s.nx,3.0))/3.0);         /*vol difference*/

 s.XP = 1.087520110412915e-01; /*Polseg vol-frac corspdng to sin.chain mix.ensemble SCF*/  
 s.Emu_s = 1.0;

 s.GL_nG=1;             /*No.of Grid points for Gauss-Laguerre quadrature*/         
 s.GL_poly = 5.0;       /*polydisperse index is (GL_poly+1.)/GL_poly */
 s.GL_Gamma = gsl_sf_gamma(s.GL_poly);     /*Gamma function using GSL library*/
 s.ct= NUMBER;  /* WEIGHT, NUMBER */
 
 /* generate the input and output file names */
 getFileName(&s);

 /* allocate memory for arrays and set initial values for the auxiliary fields */
 initialiseSystem(&s);

 /*process to solve the SCF equations */

 for(it=0;it<iter;it++)
 {
  /*To get the volume fractions from the auxiliary fields */
  getNewAvgPotential(&s); 
  getDensity(&s);
  getNewPotentialSystem(&s);

  /*calculate the field differences*/

  for(i=0;i<s.nx+1;i++)         
  {
   s.wADiff[i]=s.wNewA[i]-s.wA[i]; 
   s.wBDiff[i]=s.wNewB[i]-s.wB[i];
   s.wSDiff[i]=s.wNewS[i]-s.wS[i];
  }

  /*calculate the errors*/
  errorA = calErrors(s.wADiff,  &s);
  errorB = calErrors(s.wBDiff,  &s);      
  errorS = calErrors(s.wSDiff,  &s);
  
  freeEnergyNew=calFreeEnergy(&s);
  freeEnergyDiff=fabs(freeEnergyNew-freeEnergy);
  
  /*printing out the errors*/

  if(fmod(it,500)==0.0)
  {
   printf("%d  %e  %e  %e  %e  %.15e\n", it, errorA, errorB, errorS, freeEnergyDiff, freeEnergyNew);
  }

  if(fmod(it,5000)==0.0)
  {  
   FILE *fp;
   fp=fopen("00_iteration_1e-3","a");
   fprintf(fp,"%d  %e  %e  %e  %.15e\n", it, errorA, errorB, errorS,freeEnergyNew);
   fclose(fp);
  }
  
  if( errorA<accuracy && errorB<accuracy && errorS<accuracy ) break;

  /*set the new auxiliary potential*/
  setNewPotential(&s,it);

  freeEnergy=freeEnergyNew;

  if( fmod(it, 10000)==0 )
  {
   outputDensityEvolution(it, &s);
  }
     
  if( fmod(it, 10000)==0 )
  {  
   dataOutputField(it, &s);
  }

 } /*end of iteration*/

 getNewAvgPotential(&s);
 getDensity(&s);

  /*calculate the free energy*/
 s.freeEnergy=calFreeEnergy(&s);

  /*calculate the running time */
 finish=clock();
 runtime=(double)(finish-start)/(CLOCKS_PER_SEC);

 /*volume fractions output */
 dataOutputFile(&s, runtime);

 return 1;
}
