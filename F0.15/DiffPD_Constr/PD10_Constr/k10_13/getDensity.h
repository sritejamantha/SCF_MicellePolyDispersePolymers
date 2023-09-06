#ifndef GETDENSITYF_H
#define GETDENSITYF_H

#include "systemF.h"

double
*array1(int);  

double
**array2(int, int);

void
tridag(double *, double *, double *, double *, double *, int);

void
solvDif_CN_Neumann(double **, double *, double, double, int, int);
  
void
getDensity(systemF *);


#endif
