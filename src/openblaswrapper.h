#include "gen.h"
#include <cblas.h>


// ***********************************************************************
// Implements wrappers for various openblas functions, based on the usual
// BLAS levels. Routines included are:
//
// - a_sum:        Calculates the absolute sum of a vector. Equal to the 1-norm
// - sub_vec:      Substracts two vectors and saves it in vec1
// - copy_vec:     Copies one vector from vec to target
// - scale_vec:    Multiplies a vector with a scalar
// - dot_product:  Calculates the dotproduct between two vectors
// - dgemv_X:      Calculates the product between a mxn matrix and
//                 a m size vector. Saves the result in out. X may be
//                 'T' for transpose of the matrix or 'N' for
//                 non-transpose
// - dgemm_XY:     Calculates the product between a mxk matrix A and a 
//                 kxn matrix B and saves it in the mxn matrix C. X and Y
//                 may be 'T' or 'N' for transpose or non-transpose of A
//                 or B respectively
// - dgemm_NN_s:   Same as MatrixMult, but calcualtes A*B-C
// - matrix_trans: Calculates the transpose of an input matrix and saves it
//                 in dst
//
// - Additional Note: The library used is cblas, provided by the openblas
//                    library, which is then parallelized. Row or column/
//                    major can be chosen upon function call
//
// - Documentation about each library routine can be found on netlib.
//
// ***************************************************************************

// ----------------------------                                                         
//        BLAS LEVEl 1                                                                  
// ----------------------------                                                         
double asum(const int n, const void *vec, const char flag);
void *sub_vec(const int n, void *vec1, const void *vec2,
	      const char flag);
void *add_vec(const int n, void *vec1, const void *vec2,
	      const char flag);
void *copy_vec(const int n, const void *vec, void *target,
	       const char flag);
void *scale_vec(const int n, void *vec, double alpha,
		const char flag);
void *zscale_vec(const int n, double complex *vec,
		 double complex alpha);
double ddot_product(const int n, double *vec1, double *vec2);
/*
double zasum(const int n, double complex *vec);
void *zsub_vec(const int n, double complex *vec1, double complex *vec2);
void *zadd_vec(const int n, double complex *vec1, double complex *vec2);
void *zcopy_vec(const int n, double complex *vec, double complex *target);
void *zscale_vec(const int n, double complex *vec, double alpha);
*/
double complex zdot_product(const int n, double complex *vec1, double complex *vec2);
                                                                                        
// ----------------------------                                                         
//        BLAS LEVEL 2                                                                  
// ----------------------------                                                         
void *gemv_N(const int colA, const int rowA,
	     const void *matrix, const void *vec,
	     void *out, const char flag);
void *gemv_T(const int colA, const int rowA,
	     const void *matrix, const void *vec,
	     void *out, const char flag);
/*
void *zgemv_N(const int colA, const int rowA,
	       double complex *matrix, double complex *vec,
	       double complex *out);
void *zgemv_T(const int colA, const int rowA,
	       double complex *matrix, double complex *vec,
	       double complex *out);
*/                                                                                      
// ----------------------------                                                         
//        BLAS LEVEL 3                                                                  
// ----------------------------                                                         
void *gemm_NN(const int rowA, const int rowB, const int colA,
	      const void *matrixA, const void *matrixB,
	      void *matrixC, const char flag);         
void *gemm_NT(const int rowA, const int rowB, const int colA,
	      const void *matrixA, const void *matrixB,
	      void *matrixC, const char flag);         
void *gemm_TN(const int rowA, const int rowB, const int colA,
	      const void *matrixA, const void *matrixB,
	      void *matrixC, const char flag);         
void *gemm_TT(const int rowA, const int rowB, const int colA,
	      const void *matrixA, const void *matrixB,
	      void *matrixC, const char flag);         
void *gemm_NN_s(const int rowA, const int rowB, const int colA,
		const void *matrixA, const void *matrixB,
		void *matrixC, const char flag);       
void *gemm_NN_scal_add(const int rowA, const int colB, const int colA,
		       const void *matrixA, const void *matrixB,
		       void *matrixC, const double alpha, const double beta,
		       const char flag);
void *dmatrix_trans(const int rows, const int cols,
		    double *src, double *dst);
/*
void *zgemm_NN(const int rowA, const int rowB, const int colA,
	       double complex *matrixA, double complex *matrixB,
	       double complex *matrixC);
void *zgemm_NT(const int rowA, const int rowB, const int colA,
	       double complex *matrixA, double complex *matrixB,
	       double complex *matrixC);
void *zgemm_TN(const int rowA, const int rowB, const int colA,
	       double complex *matrixA, double complex *matrixB,
	       double complex *matrixC);
void *zgemm_NN_s(const int rowA, const int rowB, const int colA,
		 double complex *matrixA, double complex *matrixB,
		 double complex *matrixC);
*/
void *zmatrix_trans(const int rows, const int cols,
		    double complex *src, double complex *dst);

void *matrix_trans(const int rows, const int cols,
		   void *src, void *dst, const char flag);
