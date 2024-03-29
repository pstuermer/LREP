#pragma once
#include "gen.h"

// *******************************************************************
// Implements a custom COO sparse matrix storage system, mostly to
// include the Kronecker product, which is necessary to create the
// differentiation operator for higher dimensions from a one-
// dimensional setup.
//
// Includes the coo_t struct and following functions:
// - coo_malloc:      allocates memory for coo_t
// - coo_free:        frees memory for coo_t
// - coo_insert:      inserts values into a coo_t, values are required to
//                    be sorted
// - coo_kron:        calculates the kronecker product of two matrices
// - coo_get_num_row: get how many nnz values are in each row
// - coo_add:         add two coo_t matrices with different amount of
//                    nnz values
// - coo_sub:         sub two coo_t matrices with different amount of
//                    nnz values
// - coo_1norm:       calculates the 1-norm of a coo_t matrix
// *******************************************************************


typedef struct coo_t {
  int nnz, rows, cols;
  int last;
  // use only one of the two following
  double *dval;
  double complex *zval;
  
  int *col;
  int *row;
  int *numElRow;
  char flag;
} coo_t;
/*
typedef struct zcoo_t {
  int nnz, rows, cols;
  int last;
  double complex *val;
  int *col;
  int *row;
  int *numElRow;
}zcoo_t;
*/
struct coo_t *coo_malloc(const int nonz, const int numRows, const int numCols,
			 const char flag);
void coo_free(coo_t *matrix);
void coo_reset(coo_t *matrix);
void coo_insert(coo_t *matrix, double real, double imag, int colPos, int rowPos);
void dcoo_kron(coo_t *left, coo_t *right, coo_t *result);
void coo_get_num_row(coo_t *matrix);
void dcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void dcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
double coo_1norm(coo_t *matrix);

//struct zcoo_t *zcoo_malloc(const int nonz, const int numRows, const int numCols);
//void zcoo_free(zcoo_t *matrix);
//void zcoo_insert(zcoo_t *matrix, double complex val, int colPos, int rowPos);
void zcoo_kron(coo_t *left, coo_t *right, coo_t *result);
//void zcoo_get_num_row(zcoo_t *matrix);
void zcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void zcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
//double zcoo_1norm(zcoo_t *matrix);

void coo_kron(coo_t *left, coo_t *right, coo_t *result);
void coo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void coo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
