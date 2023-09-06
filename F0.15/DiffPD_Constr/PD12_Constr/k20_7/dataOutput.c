#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dataOutput.h"
#include "systemF.h"


void getFileName(systemF *s)
{
 char name[1000];
 sprintf(name,"phiP%.2f_phiS%.2f_N%d_f%.1f_xAB%.0f_xAS%.0f_xBS%.0f_kappa%.0f_dt%.3e_Lx%.0f_nx%d",s->phiPbar,s->phiSbar,s->Nn,s->f,s->xAB,s->xAS,s->xBS,s->kappa,s->dt,s->Lx,s->nx);

 strcpy(outputGeneral, "0_");
 strcat(outputGeneral, name);

 strcpy(outputDensity, "0d_");
 strcat(outputDensity, name);

 strcpy(outputField, "1_");
 strcat(outputField, name);
}


void
outputDensityEvolution(int it, systemF *s)
{
 int i;
 char c[100];
 char name[1010];
 FILE *fp;

 strcpy(name, outputDensity);

 sprintf(c,"_%d", it);
 strcat(name, c);

 fp=fopen(name, "w");
 fprintf(fp, "Title='SCMFT'\n");
 fprintf(fp, "Variables='x', 'phiA', 'phiB', 'phiS'\n");
 fprintf(fp, "ZONE T='Data'\n");
 fprintf(fp, "I=%d ZONETYPE=Ordered DATAPACKING=POINT\n", s->nx);
 for (i=0; i<s->nx+1; i++)
 {
  fprintf(fp, "%f  %.12f  %.12f  %.12f\n", s->Lx*i/s->nx, s->phiA[i], s->phiB[i], s->phiS[i]);
 }
 fclose(fp);
}



void
dataOutputFile(systemF *s, double runtime)
{
 int i;
 FILE *fp;

 fp=fopen(outputGeneral,"w");

 fprintf(fp,"the program ran %lf seconds\n", runtime);
 fprintf(fp,"the free energy of the system is %lf\n", s->freeEnergy);
 fprintf(fp,"\n");
 fprintf(fp,"\n");
 fclose(fp);


 fp=fopen(outputDensity, "w");
 fprintf(fp, "Title='SCMFT'\n");
 fprintf(fp, "Variables='x', 'phiA', 'phiB', 'phiS'\n");
 fprintf(fp, "ZONE T='Data'\n");
 fprintf(fp, "I=%d ZONETYPE=Ordered DATAPACKING=POINT\n", s->nx);

 for (i=0; i<s->nx+1; i++)
 {
  fprintf(fp, "%f  %.12f  %.12f  %.12f\n", s->Lx/s->nx*i, s->phiA[i], s->phiB[i], s->phiS[i]);
 }
 fclose(fp);


 fp=fopen(outputField,"w");
 for (i=0; i<s->nx+1; i++)
 {
  fprintf(fp, "%f  %f  %f  %f\n", s->Lx/s->nx*i, s->wA[i], s->wB[i], s->wS[i]);
 }
 fclose(fp);
}


void
dataOutputField(int it, systemF *s)
{
 int i;
 char c[100];
 char name[1010];
 FILE *fp;

 strcpy(name, outputField);

 sprintf(c,"_%d", it);
 strcat(name, c);

 fp=fopen(name, "w");
 for (i=0; i<s->nx+1; i++)
 {
   fprintf(fp, "%f  %.16f  %.16f  %.16f\n", (s->Lx/s->nx)*i, s->wA[i], s->wB[i], s->wS[i]);
 }
 fclose(fp);

}
