/* SCF in 1D Spherical coordinates ---> This is the beginning */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "systemF.h"
/*#include "randomNumbers.h"*/  /*don't know its use in this program yet */
#include "getDensity.h"
#include "energy.h"
#include "dataOutput.h"

double PI=3.14159265358979;
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

 /* set parameters for the physical system */
 /* currently all the parameters are imported from Shuanhu's work */
 /* modify them to suit to the problem of interest*/
 
 s.ns=1000;
 s.f=0.15;
 s.Nn=30;

 s.Lx=6.0*1.826;   /*radius of the sphere within which SCF is solved*/
 s.nx=80;

 s.lambda=0.01;
 s.lambda1=0.001;
 
 s.xAB=0.896;
 s.xAS=1.024*2.0;
 s.xBS=-1.0*0.05;
 s.kappa=1.171;
 s.k1 = 5.0e-03;

 s.vol=(4.0*PI*pow(s.Lx,3.0))/3.0;  /*SM -> sphere volume*/
 s.dV=s.vol/((4.0*PI*pow(s.nx,3.0))/3.0); /*SM --> vol difference*/

 s.XP = 3.0e-03;
 s.GL_nG=20;
 s.GL_poly = 5; /* only for integers, polydisperse index is (GL_poly+1.)/GL_poly */
 if (s.GL_poly==1)
 {
  s.GL_Gamma = 1.0;
 }
 else
 {
  s.GL_Gamma=1.0;
  for(i=1; i<=s.GL_poly-1; i++)
  s.GL_Gamma *= i;
 }

 s.ct= NUMBER;  /* WEIGHT, NUMBER */
 
 /* generate the input and output file names */
 getFileName(&s);

 /* allocate memory for arrays and set initial values for the auxiliary fields */
 initialiseSystem(&s);

 /*process to solve the SCF equations */

 for(it=0;it<iter;it++)
 {
  /*To get the volume fractions from the auxiliary fields */
  /*printf("%d\n \n",it);*/
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
  if(it>=0)
  {  
   errorA = calErrors(s.wADiff,  &s);
   errorB = calErrors(s.wBDiff,  &s);      
   errorS = calErrors(s.wSDiff,  &s);
  }
  else
  {
    errorA = 0.001; errorB=0.001; errorS=0.001;
  }  
  freeEnergyNew=calFreeEnergy(&s);
  freeEnergyDiff=fabs(freeEnergyNew-freeEnergy);
  
  /*printing out the errors*/

  if(fmod(it,100)==0.0)
  {
   printf("%d  %e  %e  %e  %e  %.12e\n", it, errorA, errorB, errorS, freeEnergyDiff, freeEnergyNew);
  }

  if(fmod(it,1000)==0.0)
  {  
   FILE *fp;
   fp=fopen("00_iteration_1e-3","a");
   fprintf(fp,"%d  %e  %e  %e  %e  %e %.12e\n", it, errorA, errorB, errorS, s.Qc_av,s.k1,freeEnergyNew);
   fclose(fp);
  }
  
  if( errorA<accuracy && errorB<accuracy && errorS<accuracy ) break;

  /*set the new auxiliary potential*/
  setNewPotential(&s,it);

  freeEnergy=freeEnergyNew;

  if( fmod(it, 5000)==0 )
  {
   outputDensityEvolution(it, &s);
  }
     
  if( fmod(it, 1000000)==0 )
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
