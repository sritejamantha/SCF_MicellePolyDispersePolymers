#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "energy.h"



double
calFreeEnergy( systemF *s)
{
 int i;
 double F_phi, F_kappa; /* the phi*phi term */
 double F_w; /* the phi*w term */
 double F; /* total energy */
 double coef,coef1;
 double PI=3.14159265358979;
 double F_Const=0.0;

 F_phi=0.;
 F_kappa=0.;
 F_w=0.;

 coef = 4.0*PI*pow((s->Lx/(double)s->nx),2.0); 
 coef1 = s->Lx/(double)s->nx;
 
 for (i=1; i<s->nx; i++)
 {
   F_phi += coef*pow(i,2.0)*((s->xAB*s->phiA[i]*s->phiB[i])+(s->xAS*s->phiA[i]*s->phiS[i])+(s->xBS*s->phiB[i]*s->phiS[i]));
   F_kappa += coef*pow(i,2.0)*(pow(s->phiA[i]+s->phiB[i]+s->phiS[i]-1.0, 2.0));
   F_w += coef*pow(i,2.0)*((s->phiA[i]*s->wA[i])+(s->phiB[i]*s->wB[i])+(s->phiS[i]*s->wS[i]));
 }

 F_phi += 0.5*coef*pow(s->nx,2.0)*((s->xAB*s->phiA[s->nx]*s->phiB[s->nx])+(s->xAS*s->phiA[s->nx]*s->phiS[s->nx])+(s->xBS*s->phiB[s->nx]*s->phiS[s->nx]));
 F_kappa += 0.5*coef*pow(s->nx,2.0)*(pow(s->phiA[s->nx]+s->phiB[s->nx]+s->phiS[s->nx]-1.0, 2.0));
 F_w += 0.5*coef*pow(s->nx,2.0)*((s->phiA[s->nx]*s->wA[s->nx])+(s->phiB[s->nx]*s->wB[s->nx])+(s->phiS[s->nx]*s->wS[s->nx]));
 
 F_phi *= coef1;
 F_kappa *= coef1*s->kappa/2.0;
 F_w   *= coef1;

 F_Const = 0.5*s->k1*s->vol*pow(s->phiA[10]-0.80,2.0);

 F = F_phi + F_kappa + F_Const - F_w - s->Qc_av - s->mu_s*s->Qs;

 return F;
}
