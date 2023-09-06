#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "getDensity.h"
#include "systemF.h"


double
*array1(int nx)
{
 double *m;

 m=(double *) malloc( (size_t) nx*sizeof(double) );
 return m;
}


double
**array2(int nx, int ny)
{
 int i;
 double **m;

 m=(double **) malloc( (size_t) nx*sizeof(double *) );
 m[0]=(double *) malloc( (size_t) nx*ny*sizeof(double) );

 for (i=1; i<nx; i++)
 {  
  m[i]=m[i-1]+ny;
 }
 return m;
}

void
tridag(double *a, double *b, double *c, double *r, double *U, int n)
{
 int j;
 double bet, gam[n]; 
 
 if (b[0] == 0.0) { printf("error in tridag b\n"); exit(1);}
 U[0]=r[0]/(bet=b[0]);
 for (j=1;j<n;j++)
 {
  gam[j]=c[j-1]/bet;
  bet=b[j]-a[j]*gam[j];
  if (bet == 0.0)
  {
   printf("error in tridag\n");
   exit(1);
  }
  U[j]=(r[j]-a[j]*U[j-1])/bet;
 }
 for (j=(n-2);j>=0;j--)
 {  
  U[j] -= gam[j+1]*U[j+1];
 } 
}

/*Include solvDif with Neumann bundary conditions*/

void
solvDif_CN_Neumann(double **q, double *w, double Lz, double dt, int nt, int nz)
{
 int k, it;
 double coef,coef1;
 double az[nz+1], bz[nz+1], cz[nz+1];
 double u[nz+1];

 coef=(dt*(double)(nz*nz))/(1.0*Lz*Lz);
 coef1=(dt*(double)(nz*nz))/(2.0*Lz*Lz);

 az[0]=0.0;
 bz[0]=1.0+(1.0*coef)+(dt/2.0)*w[0];
 cz[0]=-(1.0*coef);
 for (k=1; k<nz; k++)
 {
  az[k]=-(coef/2.0)+(coef1/((double)k));
  bz[k]=1.0+(coef)+(dt/2.0)*w[k];
  cz[k]=-(coef/2.0)-(coef1/((double)k));
 }
 az[nz]=-(1.0*coef);
 bz[nz]=1.0+(1.0*coef)+(dt/2.0)*w[nz];
 cz[nz]=0.0;

 u[0]=(q[0][0]*bz[0])+(q[0][1]*cz[0]);
 for(k=1; k<nz;k++)
 {  
   u[k]=(q[0][k-1]*az[k])+(q[0][k]*bz[k])+(q[0][k+1]*cz[k]);
 }
 u[nz]=(q[0][nz-1]*az[nz])+(q[0][nz]*bz[nz]);
 
 for (it=1; it<=nt; it++)
 {    
  for (k=0; k<=nz; k++)
  {  
   q[it][k]=q[it-1][k];
  } 
  for(k=0; k<=nz; k++)
  {  
    u[k]=-u[k]+(2.0*q[it][k]);
  }
  
  tridag(az, bz, cz, u, q[it], nz+1);
 }
  
 return;
}

void
getDensity(systemF *s)
{
  int i, is, iG;
  double Qc, Qs;  /* the single chain partition function */

 int nx, ns, nsA, nsB;
 double Lx;
 double f;
 double coef=0.0;
 double dt;
 int Nn,Nw=1;
 double n;
 double PI=3.141592653589793;
 Qc=0.;s->Qc_av=0.0;

 nx=s->nx;
 Lx=s->Lx;

 ns=s->ns;
 f=s->f;
 Nn=s->Nn;
 
 for(i=0;i<nx+1;i++)
 {  
  s->wShiftA[i] = s->wA[i] - s->wAAvg;
  s->wShiftB[i] = s->wB[i] - s->wBAvg;
 }
  
 for(i=0;i<nx+1;i++)
 {  
  s->phiA[i]=0.0;
  s->phiB[i]=0.0;
 }


 for(iG=0;iG<s->GL_nG; iG++)
 {  
  switch(s->ct)
  {
   case WEIGHT:
   dt = Nw*s->GL_abscissa[iG]/(s->GL_poly+1.)/ns;
   break;

   case NUMBER:
   dt = (Nn*s->GL_abscissa[iG])/(s->GL_poly*ns);
   
   break; 

   default: exit(1);
  }   

  n=0.0;
  nsA=(int)(ns*f);
  nsB=ns-nsA;

  /*Include SolvDiff functions */
  
  solvDif_CN_Neumann(s->qf, s->wShiftA, Lx, dt, nsA, nx);
  solvDif_CN_Neumann(&(s->qf[nsA]), s->wShiftB, Lx, dt, nsB, nx);
 
  solvDif_CN_Neumann(s->qd, s->wShiftB, Lx, dt, nsB, nx);
  solvDif_CN_Neumann(&(s->qd[nsB]), s->wShiftA, Lx, dt, nsA, nx);

  /*Simpson's rule for the integration*/ /*check the correctness of the implementation*/
  coef = 4.*PI*pow((s->Lx/(double)s->nx),2.0);
 
  n = (s->GL_abscissa[iG]*s->Nn)/(s->GL_poly);

  Qc=0.0;
  for(i=3;i<=nx-3;i++)
  {
   Qc+= pow(i,2.0)*s->qf[ns][i];
  }
  
  Qc+= 3.0/8.0*(pow(0.0,2.0)*s->qf[ns][0] + pow(nx,2.0)*s->qf[ns][nx]) + 7.0/6.0*(pow(1,2.0)*s->qf[ns][1]+ pow(nx-1,2.0)*s->qf[ns][nx-1]) + 23.0/24.0*(pow(2,2.0)*s->qf[ns][2] + pow(nx-2,2.0)*s->qf[ns][nx-2]);
 
  Qc *= (coef*s->Lx)/(((double)s->nx));

  s->GL_Qb[iG]=Qc;

  s->Qc_av = s->Qc_av + (s->XP*s->Emu_p[iG]*exp(-s->wAvg*n)*s->GL_Qb[iG]*s->GL_weight[iG])/(s->GL_Gamma*(double)s->Nn);


 /*calculate polymer density under potential w*/

  for (i=0; i<nx+1; i++)
  {
   s->GL_phiA[iG][i]=0.0;
   for (is=3; is<=nsA-3; is++)
   {  
    s->GL_phiA[iG][i] += s->qf[is][i]*s->qd[ns-is][i];
   }
   s->GL_phiA[iG][i] += 3./8.*(s->qf[0][i]*s->qd[ns][i] + s->qf[nsA][i]*s->qd[ns-nsA][i]) + 7./6.*(s->qf[1][i]*s->qd[ns-1][i] + s->qf[nsA-1][i]*s->qd[ns-nsA+1][i]) + 23./24.*(s->qf[2][i]*s->qd[ns-2][i] + s->qf[nsA-2][i]*s->qd[ns-nsA+2][i]);

   s->GL_phiA[iG][i] *= (s->XP*s->Emu_p[iG]*exp(-s->wAvg*n)*dt*s->GL_weight[iG])/(s->GL_Gamma*(double)s->Nn);
   
   s->phiA[i] +=s->GL_phiA[iG][i];
  }
  
  for (i=0; i<nx+1; i++)
  {
   s->GL_phiB[iG][i]=0.0;
   for (is=3; is<=nsB-3; is++)
   {  
    s->GL_phiB[iG][i] += s->qd[is][i]*s->qf[ns-is][i];
   }
   s->GL_phiB[iG][i] += 3./8.*(s->qd[0][i]*s->qf[ns][i] + s->qd[nsB][i]*s->qf[ns-nsB][i]) + 7./6.*(s->qd[1][i]*s->qf[ns-1][i] + s->qd[nsB-1][i]*s->qf[ns-nsB+1][i]) + 23./24.*(s->qd[2][i]*s->qf[ns-2][i] + s->qd[nsB-2][i]*s->qf[ns-nsB+2][i]);
   
   s->GL_phiB[iG][i] *= (s->XP*s->Emu_p[iG]*exp(-s->wAvg*n)*dt*s->GL_weight[iG])/(s->GL_Gamma*(double)s->Nn);

   s->phiB[i] +=s->GL_phiB[iG][i];
  }

 } 

 //exit(-1);
 Qs=0.;
 for (i=3; i<=nx-3; i++)
 {
   Qs += coef*pow(i,2.0)*exp(-s->wS[i]/s->Nscale);
 }
 
 Qs += 3.0/8.0*(coef*pow(0.0,2.0)*exp(-s->wS[0]/s->Nscale) + coef*pow(nx,2.0)*exp(-s->wS[nx]/s->Nscale) ) + 7.0/6.0*(coef*pow(1.0,2.0)*exp(-s->wS[1]/s->Nscale) + coef*pow(nx-1,2.0)*exp(-s->wS[nx-1]/s->Nscale) ) + 23.0/24.0*(coef*pow(2.0,2.0)*exp(-s->wS[2]/s->Nscale) + coef*pow(nx-2,2.0)*exp(-s->wS[nx-2]/s->Nscale) );

 
 Qs *= (s->Lx/(double)s->nx);
 s->Qs=Qs;

 for (i=0; i<nx+1; i++)
 {
  //s->phiS[i] = (1.0 - s->XP)*s->mu_s*exp(-s->wS[i]/s->Nscale);
  s->phiS[i] = s->Emu_s*exp(-s->wS[i]/s->Nscale);
 }

 return;
} 
