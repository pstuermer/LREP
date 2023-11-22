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
// - coo_issymmetric  checks whether a given coo_t matrix is symmetric
// - coo_kron:        calculates the kronecker product of two matrices
// - coo_get_num_row: get how many nnz values are in each row
// - coo_add:         add two coo_t matrices with different amount of
//                    nnz values
// - coo_sub:         sub two coo_t matrices with different amount of
//                    nnz values
// - coo_transpose    transposes a given coo matrix
// - coo_triu         returns the upper triangluar & diagonal part of a matrix
// - coo_1norm:       calculates the 1-norm of a coo_t matrix
// *******************************************************************


typedef struct coo_t {
  int nnz, rows, cols;    // number of nonzeros, rows and columns
  int last;               // ??
  char flag;              // signifies whether double or complex (needs rework)
  int issymetric;         // tracks whether matrix is symmetric (0=un, 1=sym)
  
  // use only one of the two following
  // either use a union or make separate coo_t structs
  double *dval;           // array of double nnz values
  double complex *zval;   // array complex double nnz values
  
  int *col;               // array of column entries
  int *row;               // array of row entries
  int *numElRow;          // ??
} coo_t;

struct coo_t *coo_malloc(const int nonz, const int numRows, const int numCols,
			 const char flag);

void coo_free(coo_t *matrix);
void coo_reset(coo_t *matrix);
void coo_insert(coo_t *matrix, double real, double imag, int colPos, int rowPos);
void coo_get_num_row(coo_t *matrix);
void coo_issymmetric(coo_t *matrix);

void dcoo_kron(coo_t *left, coo_t *right, coo_t *result);
void dcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void dcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void dcoo_triu(coo_t *matrix);

void zcoo_kron(coo_t *left, coo_t *right, coo_t *result);
void zcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void zcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void zcoo_triu(coo_t *matrix);

void coo_kron(coo_t *left, coo_t *right, coo_t *result);
void coo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void coo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR);
void coo_triu(coo_t *matrix);
double coo_1norm(coo_t *matrix);
