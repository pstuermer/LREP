#pragma once

#include "gen.h"
#include "openblaswrapper.h"
#include "physics.h"
#include "coo.h"
#include "rsbwrapper.h"
#include "splrep.h"
#include "lapackwrapper.h"

// *******************************************************************
// Implements a variety of preconditioners:
// - Conjugate Gradient
// - diagonal preconditioner
// - incomplete LU-factorization
// as well as an object that carries the choice of preconditioner and
// the preconditioner itself:
// - cond_malloc:      allocate memory for cond_t
// - cond_free:        free memory from cond_t
// - sp_setup_precond:   sets up preconditioner based on type variable
// - sp_get_diag_precond: get diagonal preconditioner from input
// - sp_conj_grad:       sparse conjugate gradient 
// - de_conj_grad:       dense conjugate gradient
// *******************************************************************


typedef struct cond_t {
  // defines the type of preconditioner to use
  int type;

  // defines the shift from 0 applied to the preconditioner
  double shift;

  // contains the diagonal inverse elements of each matrix
  // either to use as preconditioner itself or as
  // preconditioner for iterative methods
  struct rsb_mtx_t *KDiag;
  struct rsb_mtx_t *MDiag;

  // stores matrices coming from an iLU or
  // incomplete Cholesky factorization
  void *KPrecond[2];
  void *MPrecond[2];
} cond_t;
struct cond_t* cond_malloc(const int type, const double shift);

void cond_free(struct cond_t *cond);

void sp_setup_precond(struct sp_lrep_t *LREP);

void sp_get_diag_precond(struct sp_lrep_t *LREP);

void sp_conj_grad(struct rsb_mtx_t *spMatrix, double *dVec,
		  double *sol, const int size,
		  struct rsb_mtx_t *spMatrixDiag);
void sp_block_conj_gradd(struct rsb_mtx_t *spMatrix, double *dMat,
			 double *solMat, const int size, const int nrhs,
			 struct rsb_mtx_t *spMatrixDiag, const int maxIter,
			 const double tol);
void sp_block_conj_gradz(struct rsb_mtx_t *spMatrix, double complex *zMat,
			 double complex *solMat, const int size,
			 const int nrhs, struct rsb_mtx_t *spMatrixDiag,
			 const int maxIter, const double tol);

void de_conj_grad(double *deMatrix, double *dVec, double *sol,
		  const int size, struct rsb_mtx_t *deMatrixDiag);
