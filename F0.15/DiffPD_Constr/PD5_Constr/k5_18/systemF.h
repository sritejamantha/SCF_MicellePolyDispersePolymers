#ifndef SYSTEMF_H
#define SYSTEMF_H

typedef enum { UNDEFINED_INITIAL=-1, TANH, READ } initialPotential;
typedef enum { UNDEFINED_b=-1, WEIGHT, NUMBER } ChainLengthType;

typedef struct
{
  /* suppose that the monomer volumes v_0 and the kuhn lengths b for each block are the same respectively */

  int N;          /* total number of polymer segments of a diblock copolymer */
  int ns;
  int Nw,Nn;      /* weight average length and number average length */
  double phiPbar, phiSbar;

  double f;        /* segment fractions NA=N*f,  NB=N*(1-f) */
  double dt;       /* time steps in performing the EPD */
  double Dc, Ds;   /* diffusion coefficient for a single chain and a solvent molecule*/

  double lambda,lambda1;      /* parameters for convergent purpose-used in simple mixing */
  
  double Lx;       /* length of the box in x-axis scaled by Rg */
  int nx;          /* number of collocation points in the x-axis */
  
  double xAB;      /* the Flory-Huggins interaction parameter */
  double xAS;
  double xBS;
  double kappa,k1;

  double vol;         /* the volume of the box */
  double dV;
  double Qc, Qs,Qc_av;      /* single chain partition function used to calculate the free energy */
  double freeEnergy;

  initialPotential iw;
  ChainLengthType ct;
  

  double *wA, *wB, *wS;            /* the auxiliary fields */
  double *wNewA, *wNewB, *wNewS;
  double *wShiftA,*wShiftB,wAvg,wAAvg,wBAvg;
  double *wADiff, *wBDiff, *wSDiff;
  double *shift,sum;
  double *phiA, *phiB, *phiS;      /* the volume fractions */

  /* the following variables are used to solve the modified diffusion equation  */
  double **qf;
  double **qd;

  double mu_p,mu_s,phiP_Hom,phiS_Hom,*phiPA_Hom,*phiPB_Hom;
  double phiBlk_Hom, XP;
  
/********** polydispersity properties ****************/
  int GL_nG;                       /*  number of points in the Gauss-Laguerre quadrature */
  int GL_poly;                     /* polydisperse index related quantity */
  double GL_Gamma;                 /* Gamma function */
  double *GL_abscissa;             /* Gauss-Laguerre abscissas */
  double *GL_weight;               /* Gauss-Laguerre weights */
  double *GL_Qb;                   /* corresponding single chain partition function for each point */
  double **GL_phiA,**GL_phiB;
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

						 
						
