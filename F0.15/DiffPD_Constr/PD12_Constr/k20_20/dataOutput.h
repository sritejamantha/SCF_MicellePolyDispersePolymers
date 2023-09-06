#ifndef DATAOUTPUT_H
#define DATAOUTPUT_H

#include "systemF.h"

char outputGeneral[1000];  /* the volume fractions */
char outputDensity[1000];  /* the volume fractions */
char outputField[1000];   /* the fields */

void
getFileName(systemF *);

void
outputDensityEvolution(int, systemF *);

void
dataOutputFile(systemF *, double);

void
dataOutputField(int, systemF *);

#endif
