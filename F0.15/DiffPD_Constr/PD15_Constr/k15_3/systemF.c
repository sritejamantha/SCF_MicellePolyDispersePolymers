#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "systemF.h"
/*#include "randomNumbers.h"*/
#include "getDensity.h"


double 
gammln(double xx)
{
 double x,y,tmp,ser;
 static double cof[6]={76.18009172947146,-86.50532032941677, 24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
 int j;

 y=x=xx;
 tmp=x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0;j<=5;j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}


void
gaulag(double *x, double *w, int n, double alf)
{
 double EPS=3.0e-14;
 int MAXIT=1000;
 int i,its,j;
 double ai;
 double p1,p2,p3,pp,z=0.0,z1;

 for (i=0;i<n;i++)
 {
  if (i == 0)
  {
   z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
  } 
  else if (i == 1)
  {
   z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
  } 
  else
  {
   ai=i-1;
   z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
  }
  
  for (its=1;its<=MAXIT;its++)
  {
   p1=1.0;
   p2=0.0;
   for (j=1;j<=n;j++)
   {
    p3=p2;
    p2=p1;
    p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
   }
   pp=(n*p1-(n+alf)*p2)/z;
   z1=z;
   z=z1-p1/pp;
   if (fabs(z-z1) <= EPS) break;
  }
  if (its > MAXIT) { printf("too many iterations in gaulag\n"); exit(1); } 
  x[i]=z;
  w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
 }
}




void initialiseSystem(systemF *s)
{
 int i,iG;
 int ii;
 double x,phiP_in_A, phiP_in_B, phiP_in_S;
 FILE *fp1;
 fp1 = fopen("DensInp.dat","r");
 /* see getDensity.c for the details of function array3 */
 /* allocate storage space */

 
 s->wA=array1(s->nx+1);   
 s->wB=array1(s->nx+1);
 s->wS=array1(s->nx+1);
 s->wNewA=array1(s->nx+1);
 s->wNewB=array1(s->nx+1);
 s->wNewS=array1(s->nx+1);

 s->wADiff=array1(s->nx+1);
 s->wBDiff=array1(s->nx+1);
 s->wSDiff=array1(s->nx+1);

 s->wShiftA=array1(s->nx+1);
 s->wShiftB=array1(s->nx+1);
 s->shift = array1(s->GL_nG);
 
 s->wAAvg=0.;
 s->wBAvg=0.;
 s->wAvg=0.;
 s->sum=0.;
 
 s->phiPA_Hom=array1(s->nx+1);
 s->phiPB_Hom=array1(s->nx+1);
 
 s->phiA=array1(s->nx+1);
 s->phiB=array1(s->nx+1);
 s->phiS=array1(s->nx+1);
 
 for(i=0; i<s->nx+1; i++)
 {
  s->wA[i]=0.0;
  s->wB[i]=0.0;
  s->wS[i]=0.0;

  s->wNewA[i]=0.0;
  s->wNewB[i]=0.0;
  s->wNewS[i]=0.0;

  s->wADiff[i]=0.0;
  s->wBDiff[i]=0.0;
  s->wSDiff[i]=0.0;

  s->phiA[i]=0.0;
  s->phiB[i]=0.0;
  s->phiS[i]=0.0;
  
  s->wShiftA[i]=0.0;
  s->wShiftB[i]=0.0;
 }

 for(i=0;i<s->nx+1;i++)
 {
  fscanf(fp1,"%lf %lf %lf %lf",&x,&s->phiA[i],&s->phiB[i],&s->phiS[i]);
  /*s->phiA[i]=s->f*0.05;
  s->phiB[i]=(1.0 - s->f)*0.05;
  s->phiS[i]=0.95;*/
   /*printf("%d \t %.12e \t %.12e \t %.12e \t %.12e \n",i,x,s->phiA[i],s->phiB[i],s->phiS[i]);*/
  s->wA[i] = (s->xAB*s->phiB[i]) + (s->xAS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
  s->wB[i] = (s->xAB*s->phiA[i]) + (s->xBS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
  s->wS[i] = (s->xAS*s->phiA[i]) + (s->xBS*s->phiB[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
 }

 /*exit(1);*/
 fclose(fp1);



 /* note that there are NA+1 points in the contour-variable axis and nx+1 collocation points */
 s->qf=array2(s->ns+1, s->nx+1);
 s->qd=array2(s->ns+1, s->nx+1);

 for (i=0; i<s->nx+1; i++)
 {  
   s->qf[0][i]=1.;
   s->qd[0][i]=1.;   
 }

 
 s->phiP_Hom = s->phiA[s->nx] + s->phiB[s->nx];
 s->phiS_Hom = s->phiS[s->nx];
 s->phiBlk_Hom =  s->phiP_Hom +  s->phiS_Hom;

 s->mu_p = (2.*s->f*(1. - s->f)*s->xAB*s->XP) + (s->f*s->xAS + (1. - s->f)*s->xBS)*(1.0 - s->XP);

 s->mu_s = (s->f*s->xAS) + (s->xBS*(1. - s->f));
 s->mu_s = (s->mu_s)*s->XP;
 s->mu_s = (1.0 - s->XP)*exp(s->mu_s);

 /*printf("%e \t %e\t %e \t %e \t %e \t %e \t %e\n\n\n",s->vol,s->phiP_Hom,s->phiS_Hom,s->phiBlk_Hom,s->XP,s->mu_p,s->mu_s);
 exit(1);*/
 
/* caculate the abscissas and weights used in Gauss-Laguerre quadrature */
 s->GL_Qb=malloc( (size_t) s->GL_nG*sizeof(double ) );
 s->GL_abscissa=malloc( (size_t) s->GL_nG*sizeof(double ) );
 s->GL_weight=malloc( (size_t) s->GL_nG*sizeof(double ) );

 s->GL_phiA=array2(s->GL_nG, s->nx+1);
 s->GL_phiB=array2(s->GL_nG, s->nx+1);
 
 for (iG=0; iG<s->GL_nG; iG++)
 {  
  for (i=0; i<=s->nx; i++)
  {
   s->GL_phiA[iG][i]=0.0;
   s->GL_phiB[iG][i]=0.0;
  }
 } 
 gaulag( s->GL_abscissa, s->GL_weight, s->GL_nG, s->GL_poly-1.0 );

 return;
 
}

void
getNewPotentialSystem(systemF *s)
{
 int i;
 
 for(i=0;i<s->nx+1;i++)
 {
  s->k1 = 5.0e-03;
  if(i==3)
  {
   s->wNewA[i] = (s->xAB*s->phiB[i]) + (s->xAS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0)) + s->k1*s->vol*(s->phiA[i]-0.80);
   s->wNewB[i] = (s->xAB*s->phiA[i]) + (s->xBS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
   s->wNewS[i] = (s->xAS*s->phiA[i]) + (s->xBS*s->phiB[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
  }

  else
  {
   s->wNewA[i] = (s->xAB*s->phiB[i]) + (s->xAS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
   s->wNewB[i] = (s->xAB*s->phiA[i]) + (s->xBS*s->phiS[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
   s->wNewS[i] = (s->xAS*s->phiA[i]) + (s->xBS*s->phiB[i]) + (s->kappa*(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0));
  }

 }



} 

void
getNewAvgPotential(systemF *s)
{
 int i;
 double coef;
 double PI=3.14159265358979;
 int nx;
 coef = 4.*PI*pow((s->Lx/(double)s->nx),2.0);

 nx = s->nx;

 /*for(i=0;i<nx+1;i++)
 {
   printf("%.12e \t %.12e\n",s->wA[i],s->wB[i]);
 }
 exit(1);*/

 s->wAAvg=0.;
 s->wBAvg=0.;
 s->wAvg=0.;
 for(i=3;i<=nx-3;i++)
 {
  s->wAAvg += pow(i,2.0)*s->wA[i];
  s->wBAvg += pow(i,2.0)*s->wB[i];
 }

 s->wAAvg += 3.0/8.0*(pow(0.0,2.0)*s->wA[0] + pow(nx,2.0)*s->wA[nx]) + 7.0/6.0*(pow(1,2.0)*s->wA[1]+ pow(nx-1,2.0)*s->wA[nx-1]) + 23.0/24.0*(pow(2,2.0)*s->wA[2] + pow(nx-2,2.0)*s->wA[nx-2]);
 s->wAAvg = (coef*s->wAAvg/s->vol)*(s->Lx/(double)s->nx);
 
 s->wBAvg += 3.0/8.0*(pow(0.0,2.0)*s->wB[0] + pow(nx,2.0)*s->wB[nx]) + 7.0/6.0*(pow(1,2.0)*s->wB[1]+ pow(nx-1,2.0)*s->wB[nx-1]) + 23.0/24.0*(pow(2,2.0)*s->wB[2] + pow(nx-2,2.0)*s->wB[nx-2]);
 s->wBAvg = (coef*s->wBAvg/s->vol)*(s->Lx/(double)s->nx);

 s->wAvg = s->f*s->wAAvg + (1. - s->f)*s->wBAvg;

 /*printf("%.12e \t %.12e \t %.12e\n",s->wAAvg,s->wBAvg,s->wAvg);
 exit(1);*/
 return;
}

void
setNewPotential(systemF *s,int it)
{
 int i;
 double phiADiff=0.0,tol=1.0e-03;

 for (i=0; i<s->nx+1; i++)
 {
  s->wA[i] += s->lambda*s->wADiff[i];
  s->wB[i] += s->lambda*s->wBDiff[i];
  s->wS[i] += s->lambda*s->wSDiff[i];
 }

return;
}




double
calErrors(double *w,  systemF *s)
{
 int i;
 double max=0., err=0.;

 for (i=0; i<s->nx+1; i++)
 {
  err=fabs(w[i]);
  if(err>max) max=err;
 }

 return max;

}




