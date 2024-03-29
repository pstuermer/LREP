#pragma once

#include "gen.h"
#include "coo.h"
#include "grid.h"
#include "rsbwrapper.h"
#include "physics.h"


// *******************************************************************
// Implements the main linear response eigenvalue problem struct to
// be used throughout the calculation as indicated in the README.txt.
// Includes necessary routines as well:
// - splrep_malloc: allocates memory for spLREP_t struct
// - splrep_free:   frees memory of spLREP_t struct
// - setup_splrep:  sets up LREP_t struct based on input parameters
// *******************************************************************

typedef struct sp_lrep_t {
  // contains the two main matrices
  struct rsb_mtx_t *K;
  struct rsb_mtx_t *M;
  
  // contains preconditioner options
  struct cond_t *spPrecond;
  
  // eigenvectors of the current iteration
  double complex *X;
  double complex *Y;

  // eigenvectors of the previous iteration
  double complex *X1;
  double complex *Y1;

  // search directions of current iteration
  double complex *P;
  double complex *Q;

  // eigvec + search direction
  double complex *U;
  double complex *V;

  // matrix used for symmetrization
  double complex *W;

  // reduced Hamiltonian
  double complex *Hsr;

  // eigenvalues and eigenvectors
  double complex *eVecsr;
  double complex *Xsr;
  double complex *Ysr;
  double complex *eValsr;
  double complex *eVecSort;
  double complex *eValSort;

  // residual norm
  double *resNorm;

  // diagonal matrix of eigenvalues
  double complex *eigVal;

  // other stuff
  double oneNorm;
  int nev;
  int sizeSub;
  int iter;
  double tol;
  int *N;
  int size;
  int dim;
  int maxIter;
  int numComponents;
  char flag;
} sp_lrep_t;


struct sp_lrep_t *splrep_malloc(const int *N, const int nev,
				const int subspaceSize, const int dim,
				coo_t *K, coo_t *M, const int type,
				const double shift, const int maxIter,
				const double tol, const char flag);

struct sp_lrep_t *splrep_setup(char fileName[], const int *N, 
			       const double *ln, const int nev, const int subspaceSize,
			       const int dim, const double mu,
			       const int trapChoice, const double *param,
			       const int intChoice, const double *intParam,
			       const int type, const double shift,
			       const int maxIter, const double tol,
			       const int rsbTune, const char flag);

struct sp_lrep_t *splrep2C_setup(char fileName1[], char fileName2[],
				 const int *N, const double *ln,
				 const int nev, const int subspaceSize, const int dim,
				 const double mu[2], const int trapChoice,
				 const double *param, const int intChoice,
				 const double *intParam, const int type,
				 const double shift, const int maxIter,
				 const double tol, const int rsbTune,
				 const char flag);

void splrep_free(sp_lrep_t *LREP);
