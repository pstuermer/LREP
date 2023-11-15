
#include "lapackwrapper.h"

void dp_inverse(const int m, const int n, double *matrix, double *target) {
  int M, N, NRHS, LDA, LDB, RANK, LWORK, INFO;
  double *WORK, *S; 
  double RCOND, WOPT; 
  
  double *transMatrix, *transTarget;
  transMatrix = xmalloc(n*m*sizeof(double)); 
  transTarget = xmalloc(m*m*sizeof(double)); 
  matrix_trans(m, n, matrix, transMatrix, DATA_TYPE_DOUBLE); 
  matrix_trans(m, m, target, transTarget, DATA_TYPE_DOUBLE); 
  
  M = m;
  N = n;
  NRHS = M; 
  LDA = M;
  LDB = MAX(M,N); 
  RCOND = -1.0; 
  RANK = -1;
  S = xmalloc(MIN(m,n)*sizeof(double)); 
  
  // Query and allocate theoptimal workspace 
  LWORK = -1; 
  WOPT = 1.0; 
  dgelss_(&M, &N, &NRHS, transMatrix, &LDA, transTarget, 
	  &LDB, S, &RCOND, &RANK, &WOPT, &LWORK, &INFO); 
  
  // Compute Pseudoinverse and save n B 
  LWORK = WOPT; 
  WORK = xmalloc(WOPT*sizeof(double));
  dgelss_(&M, &N, &NRHS, transMatrix, &LDA, transTarget, 
	  &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO); 
  
  matrix_trans(n, m, transMatrix, matrix, DATA_TYPE_DOUBLE); 
  matrix_trans(m, m, transTarget, target, DATA_TYPE_DOUBLE); 
  
  safe_free( transTarget ); 
  safe_free( transMatrix ); 
  safe_free( WORK ); 
  
  return NULL;
}

void dcopy_matrix(const int m, const int n, 
		   double *matrix, double *target) {
  char UPLO;
  int M, N, LDA, LDB; 
  //double* awork;
  M = m;
  N = n;
  LDA = M;
  LDB = M;
  
  //awork = xmalloc(m*n*sizeof(double));
  //matrix_trans(m, n, A, awork);
  
  UPLO = 'A'; 
  
  dlacpy_(&UPLO, &M, &N, matrix, &LDA, target, &LDB);
  
  //free( awork );
  
  return NULL;
}

double dget_1_norm(const int m, const int n, double *matrix) {
  int M, N, LDA;
  double *WORK; 
  char NORM;
  double ret; 
  
  double *transMatrix;
  transMatrix = xmalloc(m*n*sizeof(double));
  dmatrix_trans(m, n, matrix, transMatrix);
  
  NORM = '1'; 
  M = m;
  N = n;
  LDA = M;
  WORK = xmalloc(M*sizeof(double)); 
  
  ret = dlange_(&NORM, &M, &N, transMatrix, &LDA, WORK);
  
  safe_free( WORK ); 
  safe_free( transMatrix );
  return ret; 
}

void dget_q(const int m, const int n, double *matrix) { 
  /* dgeqrf computes a QR factorization of a real M-by-N Matrix A */
  /**/
  /* dorgqr generates an M-by-N real matrix Q with orthonormal columns */ 
  /* which is defined as the first N columns of a product of K elementary */
  /* reflectors of order M Q = H(1)H(2)...H(k) */ 
  /* as returned by dgeqrf */ 
  
  int M, N, K, LDA, LWORK, INFO;
  double *TAU, *WORK; 
  double WKOPT; 
  //char SIDE, TRANS; 
  
  double *transMatrix;//, *transC; 
  transMatrix = xmalloc(n*m*sizeof(double)); 
  //transC = xmalloc(n*m*sizeof(double)); 
  dmatrix_trans(m, n, matrix, transMatrix); 
  //matrix_trans(m, n, C, transC); 
  
  INFO = 0; 
  M = m;
  N = n;
  LDA = M;
  WKOPT = 0.0;
  TAU = xmalloc(MIN(M,N)*sizeof(double)); 
  
  /* Query and allocate the optimal workspace */
  LWORK = -1; 
  dgeqrf_(&M, &N, transMatrix, &LDA, TAU, &WKOPT, &LWORK, &INFO);
  
  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double)); 
  
  /* Get QR-factorization */
  dgeqrf_(&M, &N, transMatrix, &LDA, TAU, WORK, &LWORK, &INFO);
  
  if (INFO < 0) { 
    printf("The algorithm failed to compute the QR-factorization.\n");
    exit( 1 );
  } 
  safe_free( WORK ); 
  
  WKOPT = 0.0;
  K = MIN(M,N); 
  INFO = 0; 
  
  // Query and allocate the optimal workspace 
  LWORK = -1; 
  dorgqr_(&M, &N, &K, transMatrix, &LDA, TAU, &WKOPT, &LWORK, &INFO);
  
  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double)); 
  
  // Get Q
  dorgqr_(&M, &N, &K, transMatrix, &LDA, TAU, WORK, &LWORK, &INFO);
  
  if (INFO < 0) { 
    printf("The algorithm failed to extract Q of a QR-factorization.\n"); 
    exit( 1 );
  } 
  
  dmatrix_trans(n, m, transMatrix, matrix); 
  // matrix_trans(n, m, transC, C);
  
  safe_free( transMatrix ); 
  //free( transC ); 
  safe_free( WORK );
  safe_free( TAU );
  
  return NULL;
}

void dcalc_sing_val(const int m, const int n, double *sVal, double *matrix) { 
  int M, N, LDA, LDU, LDVT, INFO, LWORK;
  double WKOPT; 
  double *WORK, *U, *VT;
  int *IWORK; 
  char JOBZ;
  
  double *transMatrix; 
  transMatrix = xmalloc(m*n*sizeof(double)); 
  dmatrix_trans(m, n, matrix, transMatrix); 
  
  M = m;
  N = n;
  LDA = M;
  LDU = M;
  LDVT = N; 
  INFO = 0; 
  
  JOBZ = 'N'; 
  
  U = xmalloc(LDU*M*sizeof(double));
  VT = xmalloc(LDVT*N*sizeof(double));
  IWORK = xmalloc(8*MIN(m,n)*sizeof(int));
  
  /* Query and allocate the optimal workspace */
  
  LWORK = -1; 
  dgesdd_(&JOBZ, &M, &N, transMatrix, &LDA, sVal, U, &LDU, VT, &LDVT, &WKOPT, 
	  &LWORK, IWORK, &INFO); 
  
  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double)); 
  
  /* Compute the SVD */ 
  dgesdd_(&JOBZ, &M, &N, transMatrix, &LDA, sVal, U, &LDU, VT, &LDVT, WORK, 
	  &LWORK, IWORK, &INFO); 
  
  if (INFO > 0) { 
    printf("The algorithm computing SVD failed to converge.\n");
    exit( 1 );
  } 
  
  dmatrix_trans(n, m, transMatrix, matrix); 
  
  safe_free( transMatrix ); 
  safe_free( U );
  safe_free( VT ); 
  safe_free( WORK ); 
  safe_free( IWORK );
  
  return NULL;
}

void dcalc_eig_val(const int n, double *matrix, double *eigVal, double *eigVec) {
  int N, LDA, LDVL, LDVR, INFO, LWORK;
  double WKOPT; 
  double *WORK; 
  double *WR, *WI, *VL, *VR;
  char JOBVL, JOBVR;
  
  N = n;
  LDA = N;
  LDVL = N; 
  LDVR = N; 
  WR = xmalloc(N*sizeof(double)); 
  WI = xmalloc(N*sizeof(double)); 
  VL = xmalloc(N*LDVL*sizeof(double));
  VR = xmalloc(N*LDVR*sizeof(double));
  JOBVL = 'N';
  JOBVR = 'V';
  
  double *transMatrix; 
  transMatrix = xmalloc(N*N*sizeof(double)); 
  dmatrix_trans(n, n, matrix, transMatrix); 
  
  // Query and allocate the optimal workspace 
  LWORK = -1; 
  WKOPT = 0.0;
  dgeev_(&JOBVL, &JOBVR, &N, transMatrix, &LDA, eigVal, WI, VL, &LDVL, VR, 
	 &LDVR, &WKOPT, &LWORK, &INFO); 
  
  // Calculate eigenvalues and eigenvectors 
  LWORK = (int)WKOPT; 
  WORK = xmalloc(LWORK*sizeof(double)); 
  dgeev_(&JOBVL, &JOBVR, &N, transMatrix, &LDA, eigVal, WI, VL, &LDVL, VR, 
	 &LDVR, WORK, &LWORK, &INFO); 
  /*  for(int i = 0; i < N; i++) {
    printf("%f+i*%f\t", *(eigVal+i),*(WI+i)); 
  } 
  printf("\n");
  */
  if (INFO > 0) { 
    printf("The algorithm failed to compute eigenvalues.\n"); 
    exit( 1 );
  }
  
  dmatrix_trans(n, n, transMatrix, matrix); 
  dmatrix_trans(n, n, VR, eigVec);
  
  safe_free( WR ); 
  safe_free( WI ); 
  safe_free( VL ); 
  safe_free( VR ); 
  safe_free( transMatrix );
  safe_free( WORK ); 
  return NULL;
}

void dmatrix_inverse(const int m, double *matrix) { 
  int LWORK, INFO, LDA, M; 
  double *WORK; 
  int *IPIV; 
  double WKOPT; 
  
  double *transMatrix;
  transMatrix = xmalloc(m*m*sizeof(double)); 
  dmatrix_trans(m, m, matrix, transMatrix); 
  
  M = m;
  LDA = MAX(1,M);
  INFO = 0;
  IPIV = xmalloc(M*sizeof(int));
  
  // calculate LU-factorization of matrix
  dgetrf_(&M, &M, transMatrix, &LDA, IPIV, &INFO); 
  
  // initialize workspace for inverse calculation
  LWORK = -1; 
  INFO = 0;
  WKOPT = 0.0;
  dgetri_(&M, transMatrix, &LDA, IPIV, &WKOPT, &LWORK, &INFO);
  
  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double)); 
  
  // given a LU factorization in transA, calculate the inverse
  dgetri_(&M, transMatrix, &LDA, IPIV, WORK, &LWORK, &INFO); 
  
  dmatrix_trans(M, M, transMatrix, matrix); 
  
  safe_free( transMatrix );
  safe_free( IPIV ); 
  safe_free( WORK ); 
  
  return NULL;
}


void zp_inverse(const int m, const int n, double complex *matrix, double complex *target) {
  int M, N, NRHS, LDA, LDB, RANK, LWORK, INFO;
  double complex *WORK;
  double complex *S;
  double RCOND;
  double complex WOPT, RWORK;
  
  double complex *transMatrix, *transTarget;
  transMatrix = xmalloc(n*m*sizeof(double complex));
  transTarget = xmalloc(m*m*sizeof(double complex));
  zmatrix_trans(m, n, matrix, transMatrix);
  zmatrix_trans(m, m, matrix, transTarget);
  
  M = m;
  N = n;
  NRHS = M;
  LDA = M;
  LDB = MAX(M, N);
  RCOND = -1.0;
  RANK = -1;
  S = xmalloc(MIN(m,n)*sizeof(double complex));
  
  /* Query and allocate the optimal workspace */
  LWORK = -1;
  RWORK = 5*MIN(M,N);
  WOPT = 1.0;
  zgelss_(&M, &N, &NRHS, transMatrix, &LDA, transTarget,
	  &LDB, S, &RCOND, &RANK, &WOPT, &LWORK, &RWORK, &INFO);
  
  /* Compute Pseudoinverse and save in B */
  LWORK = WOPT;
  WORK = xmalloc(WOPT*sizeof(double complex));
  zgelss_(&M, &N, &NRHS, transMatrix, &LDA, transTarget,
	  &LDB, S, &RCOND, &RANK, WORK, &LWORK, &RWORK, &INFO);
  
  zmatrix_trans(n, m, transMatrix, matrix);
  zmatrix_trans(m, m, transTarget, target);
  
  safe_free( transTarget );
  safe_free( transMatrix );
  safe_free( WORK );
  
  return NULL;
}

void zcopy_matrix(const int m, const int n, 
		   double complex *matrix, double complex *target) {
  char UPLO;
  int M, N, LDA, LDB;
  M = m;
  N = n;
  LDA = M;
  LDB = M;
  
  UPLO = 'A';
  
  zlacpy_(&UPLO, &M, &N, matrix, &LDA, target, &LDB);

  return NULL;
}

void *zget_q(const int m, const int n, double complex *matrix) {
  // zgeqrf computes a QR factorization of a real M-by-N matrix A
  //
  // zorgqr generates an M-by-N real matrix Q with orthonormal columns
  // which is defined as the first N columns of a product of K elementary
  // reflectors of order M Q = H(1)H(2)...H(k)
  // as returned by zgeqrf

  int M, N, K, LDA, LWORK, INFO;
  double complex *TAU, *WORK;
  double complex WKOPT;

  double complex *transMatrix;
  transMatrix = xmalloc(n*m*sizeof(double complex));
  zmatrix_trans(m, n, matrix, transMatrix);

  INFO = 0;
  M = m;
  N = n;
  LDA = M;
  WKOPT = 0.0;
  TAU = xmalloc(MIN(M,N)*sizeof(double complex));
  // Query and allocate the optimal workspace
  LWORK = -1;
  zgeqrf_(&M, &N, transMatrix, &LDA, TAU, &WKOPT, &LWORK, &INFO);

  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double complex));

  // get qr-factorization
  zgeqrf_(&M, &N, transMatrix, &LDA, TAU, WORK, &LWORK, &INFO);

  if(INFO < 0) {
    printf("The algorithm failed to compute the QR factorization.\n");
    exit( 1 );
  }
  safe_free( WORK );

  WKOPT = 0.0;
  K = MIN(M,N);
  INFO = 0;

  // query and allocate the optimal workspace
  LWORK = -1;
  zungqr_(&M, &N, &K, transMatrix, &LDA, TAU, &WKOPT, &LWORK, &INFO);

  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double complex));

  // get q
  zungqr_(&M, &N, &K, transMatrix, &LDA, TAU, WORK, &LWORK, &INFO);

  if (INFO < 0) {
    printf("The algorithm failed to extract Q of a QR-factorization.\n");
    exit( 1 );
  }

  zmatrix_trans(n, m, transMatrix, matrix);

  safe_free( transMatrix );
  safe_free( WORK );
  safe_free( TAU );

  return NULL;

}

void zcalc_sing_val(const int m, const int n, double *sVal, double complex *matrix) {
  int M, N, LDA, LDU, LDVT, INFO, LWORK;
  double complex WKOPT;
  double complex *WORK, *U, *VT;
  int *IWORK;
  char JOBZ;
  double *RWORK;

  double complex *transMatrix;
  transMatrix = xmalloc(m*n*sizeof(double complex));
  zmatrix_trans(m, n, matrix, transMatrix);

  M = m;
  N = n;
  LDA = M;
  LDU = M;
  LDVT = N;
  INFO = 0;

  JOBZ = 'N';
  WKOPT = 0.0;

  U = xmalloc(LDU*M*sizeof(double complex));
  VT = xmalloc(LDU*N*sizeof(double complex));
  IWORK = xmalloc(8*MIN(M,N)*sizeof(int));
  RWORK = xmalloc(5*MIN(M,N)*sizeof(double));
  printf("test\n");
  // query and allocate the optimal workspace

  LWORK = -1;
  printf("%d\n",LWORK);
  zgesdd_(&JOBZ, &M, &N, transMatrix, &LDA, sVal, U, &LDU, VT, &LDVT, &WKOPT,
	  &LWORK, RWORK, IWORK, &INFO);

  printf("test\n");
  LWORK = (int)WKOPT;
  WORK = xmalloc(LWORK*sizeof(double complex));

  // compute the svd
  zgesdd_(&JOBZ, &M, &N, transMatrix, &LDA, sVal, U, &LDU, VT, &LDVT, WORK,
	  &LWORK, RWORK, IWORK, &INFO);
  printf("test\n");

  if (INFO > 0) {
    printf("The algorithm computing a SVD failed to converge.\n");
    exit( 1 );
  }

  zmatrix_trans(n, m, transMatrix, matrix);
  safe_free( transMatrix );
  safe_free( U );
  safe_free( VT );
  safe_free( WORK );
  safe_free( IWORK );

  return NULL;

}

void zcalc_eig_val(const int n, double complex *matrix, double complex *eigVal,
		   double complex *eigVec) {
  int N, LDA, LDVL, LDVR, INFO, LWORK;
  double complex WKOPT;
  double complex *WORK;
  double complex *VL, *VR;
  char JOBVL, JOBVR;
  double *RWORK;

  N = n;
  LDA = N;
  LDVL = N;
  LDVR = N;
  VL = xmalloc(N*LDVL*sizeof(double complex));
  VR = xmalloc(N*LDVR*sizeof(double complex));
  RWORK = xmalloc(2*N*sizeof(double));
  JOBVL = 'N';
  JOBVR = 'V';

  double complex *transMatrix;
  transMatrix = xmalloc(N*N*sizeof(double complex));
  zmatrix_trans(n, n, matrix, transMatrix);

  // querz and allocate the optimal workspace
  LWORK = -1;
  WKOPT = 0.0;
  zgeev_(&JOBVL, &JOBVR, &N, transMatrix, &LDA, eigVal, VL, &LDVL, VR,
	 &LDVR, &WKOPT, &LWORK, RWORK, &INFO);

  // calculate eigenvalues and eigenvectors
  LWORK = (int)WKOPT;
  WORK = xmalloc(LWORK*sizeof(double complex));
  zgeev_(&JOBVL, &JOBVR, &N, transMatrix, &LDA, eigVal, VL, &LDVL, VR,
	 &LDVR, WORK, &LWORK, RWORK, &INFO);

  /*  for(int i = 0; i < n; i++) {
    printf("%.7e %.7e\t",creal(eigVal[i]),cimag(eigVal[i]));
  }
  printf("\n");
  */
  if (INFO>0) {
    printf("The algorithm failed to compute eigenvalues.\n");
    exit( 1 );
  }
  zmatrix_trans(n, n, transMatrix, matrix);
  zmatrix_trans(n, n, VR, eigVec);

  safe_free( VL );
  safe_free( VR );
  safe_free( transMatrix );
  safe_free( WORK );
  safe_free( RWORK );
  return NULL;
}

void zmatrix_inverse(const int m, double complex *matrix) {
  int LWORK, INFO, LDA, M;
  double complex *WORK;
  int *IPIV;
  double complex WKOPT;

  double complex *transMatrix;
  transMatrix = xmalloc(m*m*sizeof(double complex));
  zmatrix_trans(m, m, matrix, transMatrix);

  M = m;
  LDA = MAX(1,M);
  INFO = 0;
  IPIV = xmalloc(M*sizeof(int));

  // calculate LU-factorization of matrix
  zgetrf_(&M, &M, transMatrix, &LDA, IPIV, &INFO);

  // initialize workspace for inverse calculation

  LWORK = -1;
  INFO = 0;
  WKOPT = 0.0;
  zgetri_(&M, transMatrix, &LDA, IPIV, &WKOPT, &LWORK, &INFO);

  LWORK = WKOPT;
  WORK = xmalloc(LWORK*sizeof(double complex));

  // given a LU-factorization in transMatrix, calculate the inverse
  zgetri_(&M, transMatrix, &LDA, IPIV, WORK, &LWORK, &INFO);

  zmatrix_trans(m, m, transMatrix, matrix);

  safe_free( transMatrix );
  safe_free( IPIV );
  safe_free( WORK );

  return NULL;
}



void p_inverse(const int m, const int n, void *matrix,
		void *target, const char flag) {

  if(flag == 'D') {
    dp_inverse(m, n, matrix, target);
  } else if (flag == 'Z') {
    zp_inverse(m, n, matrix, target);
  }

  return NULL;
}

void copy_matrix(const int m, const int n, void *matrix,
		  void *target, const char flag) {

  if(flag == 'D') {
    dcopy_matrix(m, n, matrix, target);
  } else if(flag == 'Z') {
    zcopy_matrix(m, n, matrix, target);
  }

  return NULL;
}

void get_q(const int m, const int n, void *matrix,
	    const char flag) {
  if(flag == 'D') {
    dget_q(m, n, matrix);
  } else if(flag == 'Z') {
    zget_q(m, n, matrix);
  }

  return NULL;
}

void calc_sing_val(const int m, const int n, double *sVal, void *matrix,
		    const char flag) {
  if(flag == 'D') {
    dcalc_sing_val(m, n, sVal, matrix);
  } else if(flag =='Z') {
    zcalc_sing_val(m, n, sVal, matrix);
  }

  return NULL;
}

void calc_eig_val(const int m, void *matrix, void *eigVal, void *eigVec,
		   const char flag) {
  if(flag == 'D') {
    dcalc_eig_val(m, matrix, eigVal, eigVec);
  } else if(flag == 'Z') {
    zcalc_eig_val(m, matrix, eigVal, eigVec);
  }

  return NULL;
}

void matrix_inverse(const int m, void *matrix, const char flag) {
  if (flag == 'D') {
    dmatrix_inverse(m, matrix);
  } else if(flag == 'Z') {
    zmatrix_inverse(m, matrix);
  }

  return NULL;
}

double z_norm(const int m, const int n, double complex *matrix, char NORM) {                      
  int M, N, LDA;                                                                                  
  double *WORK;                                                                                   
  double ret;                                                                                     
                                                                                                  
  double complex *transMatrix;                                                                    
  transMatrix = xmalloc(n*m*sizeof(double complex));                                              
  zmatrix_trans(m, n, matrix, transMatrix);                                                       
                                                                                                  
  M = m;                                                                                          
  N = n;                                                                                          
  LDA = M;                                                                                        
  WORK = xmalloc(M*sizeof(double));                                                               
                                                                                                  
  ret = zlange_(&NORM, &M, &N, transMatrix, &LDA, WORK);                                          
                                                                                                  
  safe_free( WORK );                                                                              
  safe_free( transMatrix );                                                                       
                                                                                                  
  return ret;                                                                                     
}                                                                                                 
                                                                                                  
double d_norm(const int m, const int n, double *matrix, char NORM) {                              
  int M, N, LDA;                                                                                  
  double *WORK;                                                                                   
  double ret;                                                                                     
                                                                                                  
  double *transMatrix;                                                                            
  transMatrix = xmalloc(n*m*sizeof(double));                                                      
  dmatrix_trans(m, n, matrix, transMatrix);                                                       
                                                                                                  
  M = m;                                                                                          
  N = n;                                                                                          
  LDA = M;                                                                                        
  WORK = xmalloc(M*sizeof(double));                                                               
                                                                                                  
  ret = dlange_(&NORM, &M, &N, transMatrix, &LDA, WORK);                                          
                                                                                                  
  safe_free( WORK );                                                                              
  safe_free( transMatrix );                                                                       
                                                                                                  
  return ret;                                                                                     
}                                                                                                 
                                                                                                  
double matrix_norm(const int m, const int n, void *matrix, const char NORM,                       
		   const char flag) {                                                              
  double ret;                                                                                     
  ret = -1.0;                                                                                     
                                                                                                  
  if(flag == 'D') {                                                                               
    ret = d_norm(m, n, matrix, NORM);                                                             
  } else if(flag == 'Z') {                                                                        
    ret = z_norm(m, n, matrix, NORM);                                                             
  }                                                                                               
                                                                                                  
  return ret;                                                                                     
}
  
