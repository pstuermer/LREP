
#include "openblaswrapper.h"

// ----------------------------
//  BLAS Level 1
// ----------------------------
 
double asum(const int n, const void *vec, const char flag) {
  double ret = 0.0;
  int INCX;
  INCX = 1;
 
  if(flag == 'D') {
    ret = cblas_dasum(n, vec, INCX);
  } else if (flag == 'Z') {
    ret = cblas_dzasum(n, vec, INCX);
  }
  return ret;
}
 
void *sub_vec(const int n, void *vec1, const void *vec2,
	      const char flag) { 
  double DA[2];
 
  int INCX, INCY;  
  INCX = 1;  
  INCY = 1;  

  DA[0] = -1.0; 
  DA[1] = 0.0;

  if(flag == 'D') {
    cblas_daxpy(n, DA[0], vec2, INCX, vec1, INCY);
  } else if (flag == 'Z') {
    cblas_zaxpy(n, DA, vec2, INCX, vec1, INCY);
  }
 
  return NULL;  
}
 
void *add_vec(const int n, void *vec1, const void *vec2,
	      const char flag) {  
  double DA[2]; 

  int INCX, INCY;  
  INCX = 1;  
  INCY = 1;  

  DA[0] = 1.0;  
  DA[1] = 0.0;

  if(flag == 'D') {
    cblas_daxpy(n, DA[0], vec2, INCX, vec1, INCY);
  } else if (flag == 'Z') {
    cblas_zaxpy(n, DA, vec2, INCX, vec1, INCY);
  }
 
  return NULL;  
}
 
void *copy_vec(const int n, const void *vec, void *target,
	       const char flag) {
  int INCX, INCY;  
  INCX = 1;  
  INCY = 1;  

  if(flag == 'D') {
    cblas_dcopy(n, vec, INCX, target, INCY); 
  } else if(flag == 'Z') {
    cblas_zcopy(n, vec, INCX, target, INCY);
  }
 
  return NULL;  
}
 
void *scale_vec(const int n, void *vec, double alpha,
		const char flag) { 
  int INCX;  
  INCX = 1;  

  if(flag == 'D') {
    cblas_dscal(n, alpha, vec, INCX);  
  } else if(flag == 'Z') {
    cblas_zdscal(n, alpha, vec, INCX);
  }
 
  return NULL;  
}

void *zscale_vec(const int n, double complex *vec,
		 double complex alpha) {
  int INCX;
  INCX = 1;
  

  cblas_zscal(n, &alpha, vec, INCX);

  return NULL;
}
 
double ddot_product(const int n, double *vec1, double *vec2) {
  double ret = 0.0;
  int INCX, INCY;  
  INCX = 1;  
  INCY = 1;  
 
  ret = cblas_ddot(n, vec1, INCX, vec2, INCY);
  return ret;
}

/*
double zasum(const int n, double complex *vec) {
  double complex ret = 0.0;
  int INCX;
  INCX = 1;
  ret = cblas_dzasum(n, vec, INCX);
  return ret;
}

void *zsub_vec(const int n, double complex *vec1, double complex *vec2) {
  double DA[2];
  int INCX, INCY;
  INCX = 1;
  INCY = 1;
  DA[0] = -1.0;
  DA[1] = 0.0;

  cblas_zaxpy(n, DA, vec2, INCX, vec1, INCY);

  return NULL;
}

void *zadd_vec(const int n, double complex *vec1, double complex *vec2) {
  double DA[2];
  int INCX, INCY;
  INCX = 1;
  INCY = 1;
  DA[0] = 1.0;
  DA[1] = 0.0;

  cblas_zaxpy(n, DA, vec2, INCX, vec1, INCY);

  return NULL;
}

void *zcopy_vec(const int n, double complex *vec, double complex *target) {
  int INCX, INCY;
  INCX = 1;
  INCY = 1;
  cblas_zcopy(n, vec, INCX, target, INCY);

  return NULL;
}

void *zscale_vec(const int n, double complex *vec, double alpha) {
  int INCX;
  INCX = 1;
  cblas_zscal(n, &alpha, vec, INCX);

  return NULL;
}
*/
double complex zdot_product(const int n, double complex *vec1, double complex *vec2) {
  double complex ret = 0.0;
  int INCX, INCY;
  INCX = 1;
  INCY = 1;

  cblas_zdotc_sub(n, vec1, INCX, vec2, INCY, &ret);

  return ret;
}

// ----------------------------
//  BLAS Level 2
// ----------------------------
 
void *gemv_N(const int colA, const int rowA, const void *matrix, 
	     const void *vec, void *out, const char flag) {  
 
  int lda, incx, incy;
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, colA); 
  incx = 1;  
  incy = 1;  

  if(flag == 'D') {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rowA, colA, one[0],
		matrix, lda, vec, incx, zero[0], out, incy);
  } else if(flag == 'Z') {
    cblas_zgemv(CblasRowMajor, CblasNoTrans, rowA, colA, one,
		matrix, lda, vec, incx, zero, out, incy);
  }

  return NULL;  
}
 
void *gemv_T(const int colA, const int rowA, const void *matrix,
	     const void *vec, void *out, const char flag) {  
  int lda, incx, incy;
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, colA); 
  incx = 1;  
  incy = 1;  
 
  if(flag == 'D') {
    cblas_dgemv(CblasRowMajor, CblasTrans, rowA, colA, one[0], 
		matrix, lda, vec, incx, zero[0], out, incy); 
  } else if(flag == 'Z') {
    cblas_zgemv(CblasRowMajor, CblasTrans, rowA, colA, one,
		matrix, lda, vec, incx, zero, out, incy);
  }

  return NULL;  
}
/*
void *zgemv_N(const int colA, const int rowA,
	 double complex *matrix, double complex *vec,
	 double complex *out) {

  int lda, incx, incy;
  double one[2], zero[2];
  lda = MAX(1, colA);
  incx = 1;
  incy = 1;
  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[0] = 0.0;

  cblas_zgemv(CblasRowMajor, CblasNoTrans, rowA, colA, one,
	      matrix, lda, vec, incx, zero, out, incy);
  return NULL;
}

void *zgemv_T(const int colA, const int rowA,
	double complex *matrix, double complex *vec,
	double complex *out) {

  int lda, incx, incy;
  double one[2], zero[2];
  lda = MAX(1, colA);
  incx = 1;
  incy = 1;
  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;
  
  cblas_zgemv(CblasRowMajor, CblasTrans, rowA, colA, one,
	      matrix, lda, vec, incx, zero, out, incy);

  return NULL;
}
*/
// ----------------------------
//  BLAS Level 3
// ----------------------------
 
void *gemm_NN(const int rowA, const int rowB, const int colA,
	      const void *matrixA, const void *matrixB, 
	      void *matrixC, const char flag) {  

  int lda, ldb, ldc;  
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, colA); 
  ldb = MAX(1, rowB); 
  ldc = MAX(1, rowB); 

  if(flag == 'D') {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, rowB, colA,  
		one[0], matrixA, lda, matrixB, ldb, zero[0], matrixC, ldc); 
  } else if(flag == 'Z') {
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, rowB, colA,
		one, matrixA, lda, matrixB, ldb, zero, matrixC, ldc);
  }

  return NULL;  
}
 
void *gemm_NT(const int rowA, const int rowB, const int colA, 
	      const void *matrixA, const void *matrixB, 
	      void* matrixC, const char flag) {  
  int lda, ldb, ldc;  
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, colA); 
  ldb = MAX(1, colA); 
  ldc = MAX(1, rowB); 

  if(flag == 'D') {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowA, rowB, colA, 
		one[0], matrixA, lda, matrixB, ldb, zero[0], matrixC, ldc); 
  } else if(flag == 'Z') {
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rowA, rowB, colA,
		one, matrixA, lda, matrixB, ldb, zero, matrixC, ldc);
  }

  return NULL;  
}
 
void *gemm_TN(const int rowA, const int colB, const int colA,
	       const void *matrixA, const void *matrixB, 
	       void *matrixC, const char flag) {  
  int lda, ldb, ldc;  
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, rowA); 
  ldb = MAX(1, colB); 
  ldc = MAX(1, colB); 

  if(flag == 'D') {
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, rowA, colB, colA, 
		one[0], matrixA, lda, matrixB, ldb, zero[0], matrixC, ldc);
  } else if(flag == 'Z') {
    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, rowA, colB, colA,
		one, matrixA, lda, matrixB, ldb, zero, matrixC, ldc);
  }

  return NULL;  
}

void *gemm_TT(const int rowA, const int rowB, const int colA,
	       const void *matrixA, const void *matrixB, 
	       void *matrixC, const char flag) {  
  int lda, ldb, ldc;  
  double one[2], zero[2];

  one[0] = 1.0;
  one[1] = 0.0;
  zero[0] = 0.0;
  zero[1] = 0.0;

  lda = MAX(1, rowA); 
  ldb = MAX(1, colA); 
  ldc = MAX(1, rowB); 

  if(flag == 'D') {
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, rowA, rowB, colA,
		one[0], matrixA, lda, matrixB, ldb, zero[0], matrixC, ldc);  
  } else if(flag == 'Z') {
    cblas_zgemm(CblasRowMajor, CblasTrans, CblasTrans, rowA, rowB, colA,
		one, matrixA, lda, matrixB, ldb, zero, matrixC, ldc);
  }

  return NULL;  
}
 
void *gemm_NN_s(const int rowA, const int rowB, const int colA, 
		const void *matrixA, const void *matrixB,
		void *matrixC, const char flag) {
  int lda, ldb, ldc;  
  double one[2], minus[2];

  one[0] = 1.0;
  one[1] = 0.0;
  minus[0] = -1.0;
  minus[1] = 0.0;

  lda = MAX(1, colA); 
  ldb = MAX(1, rowB); 
  ldc = MAX(1, rowB); 
 
  if(flag == 'D') {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, rowB, colA,  
		one[0], matrixA, lda, matrixB, ldb, minus[0], matrixC, ldc);
  } else if(flag == 'Z') {
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, rowB, colA,
		one, matrixA, lda, matrixB, ldb, minus, matrixC, ldc);
  }
  return NULL;  
}

void *gemm_NN_scal_add(const int rowA, const int colB, const int colA,                            
                       const void *matrixA, const void *matrixB,                                  
                       void *matrixC, const double alpha, const double beta,                      
                       const char flag) {                                                         
                                                                                                  
  int lda, ldb, ldc;                                                                              
  double dalpha[2], dbeta[2];                                                                     
                                                                                                  
  dalpha[0] = alpha;                                                                              
  dalpha[1] = 0.0;                                                                                
  dbeta[0] = beta;                                                                                
  dbeta[1] = 0.0;                                                                                 
                                                                                                  
  lda = MAX(1, colA);                                                                             
  ldb = MAX(1, colB);                                                                             
  ldc = MAX(1, colB);                                                                             
                                                                                                  
  if(flag == 'D') {                                                                               
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA,                      
                dalpha[0], matrixA, lda, matrixB, ldb, dbeta[0], matrixC, ldc);                   
  } else if(flag == 'Z') {                                                                        
    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA,                      
                dalpha, matrixA, lda, matrixB, ldb, dbeta, matrixC, ldc);                         
  }                                                                                               
  return NULL;                                                                                    
}
 
void *dmatrix_trans(const int rows, const int cols, 
		    double* src, double* dst) { 
  // pragma omp parallel for
  // rows and cols correspond to src, not dst 
  int i, j;  
  for (int k = 0; k < rows*cols; k++) { 
    i = k/rows; 
    j = k%rows; 
    dst[k] = src[cols*j+i]; 
  } 
 
  return NULL;  
}

void *zmatrix_trans(const int rows, const int cols,
		    double complex *src, double complex *dst) {
  // rows and cols correspond to src, not dst
  int i, j;
  for(int k = 0; k < rows*cols; k++) {
    i = k/rows;
    j = k%rows;
    dst[k] = src[cols*j+i];
  }

  return NULL;
}

void *matrix_trans(const int rows, const int cols,
		   void *src, void *dst, const char flag) {
  if(flag == 'D') {
    printf("picked double transpose.\n");
    dmatrix_trans(rows, cols, src, dst);
  } else if (flag == 'Z') {
    zmatrix_trans(rows, cols, src, dst);
  }

  return NULL;
}
