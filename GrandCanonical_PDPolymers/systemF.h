#ifndef SYSTEMF_H
#define SYSTEMF_H

typedef enum { UNDEFINED_INITIAL=-1, TANH, READ } initialPotential;
typedef enum { UNDEFINED_b=-1, WEIGHT, NUMBER } ChainLengthType;

typedef struct
{
  /* suppose that the monomer volumes v_0 and the kuhn lengths b for each block are the same respectively */

  int ns;          /*Chain grid points to solve diffusion equation*/
  int Nn;          /* Scaled Number average length. set to 1 */
  double Nscale;   /*Chain length scaling factor*/

  double f;        /* Fraction of solvophobic groups in a chain */

  double lambda;   /* parameter to update potential-used in simple mixing */
  
  double Lx;       /* Radius of the spherical box [scaled by Rg(Nscale.Nn)] */
  int nx;          /* number of spatial grid points*/
  
  double xAB;      /* the Flory-Huggins interaction parameter */
  double xAS;
  double xBS;
  double kappa;    /*Parameter to control compressibility*/

  double vol,dV;
  double Qc,Qs,Qc_av; /* partition function used to calculate the free energy */
  double freeEnergy;

  initialPotential iw;
  ChainLengthType ct;
  

  double *wA, *wB, *wS;            /* the auxiliary fields */
  double *wNewA, *wNewB, *wNewS;
  double *wShiftA,*wShiftB,wAvg,wAAvg,wBAvg;
  double *wADiff, *wBDiff, *wSDiff;
  double *phiA, *phiB, *phiS;      /* the volume fractions */

  /* the following variables are used to solve the modified diffusion equation  */
  double **qf;
  double **qd;

  double *Emu_p,Emu_s;    /*Exponential of polymer and solvent chem pot*/
  double XP;      /*Pol.Seg Vol.Frac used to compute Emu_p from Sing.Chn Mixed.Ensemble SCF*/
  
/********** polydispersity properties ****************/
  int GL_nG;                       /*  number of points in the Gauss-Laguerre quadrature */
  double GL_poly;                  /* polydisperse index related quantity */
  double GL_Gamma;                 /* Gamma function */
  double *GL_abscissa;             /* Gauss-Laguerre abscissas */
  double *GL_weight;               /* Gauss-Laguerre weights */
  double *GL_Qb;                   /* corresponding single chain partition function for each point */
  double **GL_phiA,**GL_phiB;      /*Polymer segment volume fractions*/
}systemF;


double 
gammln(double);

void
gaulag(double *, double *, int, double);

void
initialiseSystem(systemF *);

double
calErrors(double *,  systemF *);

void
getNewPotentialSystem(systemF *);

void
getNewAvgPotential(systemF *);

void
setNewPotential(systemF *,int);

#endif

						 
						
