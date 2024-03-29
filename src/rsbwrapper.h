#pragma once
#include "gen.h"
#include "coo.h"

// *******************************************************************
// Includes wrapper functions for the Recursive-sparse-block matrix 
// library, which is used for sparse-dense matrix-matrix and sparse-dense
// matrix-vector multiplication.
//
// Routines include:
// - rsb_init:         initializes the library using some tuning parameter
// - rsb_mtx_from_coo: creates a rsb matrix from the custom coo_t format
// - rsb_SPMM:         sparse-dense matrix-matrix multiplication
// - rsb_SPMM_sub      sparse-dense matrix-matrix multiplication with
//                     C = C - A*B
// - rsb_SPMV:         sparse-dense matrix-vector multiplication
// - rsb_SPVM:         dense-sparse vector-matrix multiplication
// - rsb_SPSM:         sparse-dense matrix-matrix multiplication with
//                     A upper or lower triangular
// - rsb_SPSV:         sparse-dense matrix-vector multiplication with
//                     A upper or lower triangular
// - rsb_get_prec:     calculates incomplete LU-factorization from input
// - rsb_tune_SPMM:    optimizes thread count and recursive block
//                     structure to speed-up matrix-matrix multiplication.
//                     Considering SPMM takes up to 80% of the runtime,
//                     this is absolutely crucial
// *******************************************************************

void rsb_init(const int rsbTune);
struct rsb_mtx_t *rsb_mtx_from_coo(coo_t *matrix);
void rsb_SPMM(struct rsb_mtx_t *spMatrix, const void *dMatrix,
	      void *res, const int nev, const char flag);
void rsb_SPMM_sub(struct rsb_mtx_t *spMatrix, const void *dMatrix,
		  void *res, const int nev, const char flag);
void rsb_SPMM_scal_add(struct rsb_mtx_t *spMatrix, const void *dMatrix,
		       void *res, int nev, const double alpha,
		       const double beta, const char flag);
void rsb_SPMV(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag);
void rsb_SPVM(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag);
void rsb_SPSM(struct rsb_mtx_t *spMatrix, const void *dMatrix,
	      void *res, const int nev, const char flag);
void rsb_SPSV(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag);
void rsb_get_prec(struct rsb_mtx_t *spMatrix, void *opdp[2]);
void* rsb_tune_SPMM(struct rsb_mtx_t *spMatrix, const double *dMatrix, 
		   const int nev, const int size, const int tn,
		   const char flag);
/*
struct rsb_mtx_t *rsb_mtx_from_zcoo(zcoo_t *matrix);
void rsb_SPMM_test(struct rsb_mtx_t *spMatrix, const void *zMatrix,
		   double complex *res, const int nev);
*/
