
#include "precond.h"
#include <rsb.h>

struct cond_t *cond_malloc(const int type, const double shift) {
  cond_t *cond = xmalloc(sizeof(struct cond_t));

  cond -> type = type;
  cond -> shift = shift;

  return cond;
}

void cond_free(struct cond_t *cond) {
  // check if rsb matrices exist or not (need to chack that in the docs)
  const int type = cond->type;

  if(!type) {
    // do nothing, as no preconditioner chosen
    ;
  } else if (type == 1) {
    rsb_mtx_free( cond->MDiag );
    rsb_mtx_free( cond->KDiag );
  } else if (type == 2) {
    rsb_mtx_free( cond->MDiag );
    rsb_mtx_free( cond->KDiag );
  } else if (type == 3) {
    rsb_mtx_free( cond->MPrecond[0] );
    rsb_mtx_free( cond->MPrecond[1] );
    rsb_mtx_free( cond->KPrecond[0] );
    rsb_mtx_free( cond->KPrecond[1] );
  } else {
    fprintf(stderr, "Preconditioner type not recognized.\n");
    safe_free( cond );
    exit(1);
  }
 
  safe_free( cond );
}

void* sp_setup_precond(struct sp_lrep_t *LREP) {
  printf("get precond.\n");
  // what do I want to do here?  
  // at end of setup of my LREP setup, I go through the chosen options   
  // of the preconditioner (which currently only allows ConjGrad anyways 
  // and set that one up. Later I need something similar, but for   
  // applyPrecond. (Potentially could also use the diag precond)    
   
  // either switch or if    
  int type = LREP->spPrecond->type;
   
  if (!type) {    
    // do nithing, as no preconditioner chosen  
    ;
  } else if (type == 1) {   
    // diagonal preconditioner   
    sp_get_diag_precond(LREP);
  } else if (type == 2) {   
    // conjugate gradient preconditioner
    sp_get_diag_precond(LREP);
  } else if (type == 3) {   
    // ILU(0) preconditioner (doesnt have shift in it yet)
    // also generally not working
    rsb_get_prec(LREP->K, LREP->spPrecond->KPrecond);
    rsb_get_prec(LREP->M, LREP->spPrecond->MPrecond);
  } else {   
    fprintf(stderr, "Preconditioner type not recognized.\n");
    exit(1);
  }
  
  return NULL;
}

void* sp_get_diag_precond(struct sp_lrep_t *LREP) { 
  // do it for both K and M 
  int size = LREP->size;
  coo_t *MDiag, *KDiag;
  MDiag = coo_malloc(size, size, size, LREP->flag);
  KDiag = coo_malloc(size, size, size, LREP->flag);

  double complex *M, *K;
  M = xmalloc(size*sizeof(double complex));
  K = xmalloc(size*sizeof(double complex));
  
  // K first 
  rsb_mtx_get_vec(LREP->K, K, RSB_EXTF_DIAG);
  for(int i = 0; i < size; i++)
    coo_insert(KDiag, 1.0/(creal(K[i]-LREP->spPrecond->shift)),
	       0.0, i, i);
  LREP->spPrecond->KDiag = rsb_mtx_from_coo(KDiag);


  // M next  
  rsb_mtx_get_vec(LREP->M, M, RSB_EXTF_DIAG);
  for (int i = 0; i < size; i++) 
    coo_insert(MDiag, 1.0/(creal(M[i]-LREP->spPrecond->shift)),
	       0.0, i, i);
  LREP->spPrecond->MDiag = rsb_mtx_from_coo(MDiag);

  
  coo_free( KDiag );
  safe_free( K );
  
  coo_free( MDiag );
  safe_free( M );
  
  return NULL;
}

void* de_conj_grad(double *deMatrix, double *dVec, double *sol, 
		   const int size, struct rsb_mtx_t *deMatrixDiag) {   
  // solves the symmetric positive definite linear system Ax=b 
  // using the Conjugate Gradient method (with preconditioning)
  int maxIter;
  double tol;
  double err;
  
  double alpha, rho, rho1, beta;
  const char flag = DATA_TYPE_DOUBLE;
  
  maxIter = 100;
  tol = 10e-8;
  
  // get norm of b
  // if norm of b == 0, return null vector 
  double bNorm2 = sqrt(ddot_product(size, dVec, dVec));
  
  if (bNorm2 == 0.0) { 
    bNorm2 = 1.0;
    err = 0.0;
    memset(sol, 0, size*sizeof(double));
    
    return NULL;
  }
  
  double *r, *z, *p;
  r = xmalloc(size*sizeof(double));
  z = xmalloc(size*sizeof(double));
  p = xmalloc(size*sizeof(double));
  copy_vec(size, dVec, r, flag);
  
  double *work;
  work = xmalloc(size*sizeof(double));
  
  // calculate residual r   
  gemv_T(size, size, deMatrix, sol, work, flag);
  
  sub_vec(size, r, work, flag);
  err = sqrt(ddot_product(size, r, r))/bNorm2;
  
  if (err <= tol) 
    return NULL;

  for (int i = 0;i < maxIter;i++) { 
    rsb_SPMV(deMatrixDiag, r, z, flag);
    rho = ddot_product(size, r, z);
   
    if (i) { 
      beta = rho/rho1;
      scale_vec(size, p, beta, flag);
      add_vec(size, p, z, flag);
    } else { 
      copy_vec(size, z, p, flag);
    }   
    
    gemv_T(size, size, deMatrix, p, work, flag);
    alpha = rho/ddot_product(size, p, work);
    scale_vec(size, p, alpha, flag);
    add_vec(size, sol, p, flag);
    
    scale_vec(size, work, alpha, flag);
    sub_vec(size, r, work, flag);
    err = sqrt(ddot_product(size, r, r))/bNorm2;
    
    if (err <= tol)    
      return NULL;
    
    
    rho1 = rho;
  }
  
  safe_free( r );
  safe_free( z );
  safe_free( p );
  
  return NULL;
}

void *sp_conj_grad(struct rsb_mtx_t *spMatrix, double *dVec,
		   double *sol, const int size, 
		   struct rsb_mtx_t *spMatrixDiag) { 
  // solves the symmetric positive definite linear system Ax=b 
  // using the Conjugate Gradient method (with preconditioning)
  int maxIter;
  double tol;
  double err;
  
  double alpha, rho, rho1, beta;
  const char flag = DATA_TYPE_DOUBLE;
  
  maxIter = 100;
  tol = 10e-8;
  
  // get norm of b
  // if norm of b == 0, return null vector
  double bNorm2;
  bNorm2 = sqrt(ddot_product(size, dVec, dVec));
  
  if (bNorm2 == 0.0) { 
    bNorm2 = 1.0;
    err = 0.0;
    memset(sol, 0, size*sizeof(double));
    
    return NULL;
  }
  
  double *r, *z, *p;
  r = xmalloc(size*sizeof(double));
  z = xmalloc(size*sizeof(double));
  p = xmalloc(size*sizeof(double));
  copy_vec(size, dVec, r, flag);
  
  double *work;
  work = xmalloc(size*sizeof(double));
  
  // calculate residual r
  rsb_SPMV(spMatrix, sol, work, flag);
  sub_vec(size, r, work, flag);
  
  err = sqrt(ddot_product(size, r, r))/bNorm2;

  
  if (err <= tol) 
    return NULL;
  for(int i = 0;i < maxIter;i++) {
    rsb_SPVM(spMatrixDiag, r, z, flag);
    rho = ddot_product(size, r, z);

    if( i > 1 ) { 
      beta = rho/rho1;
      scale_vec(size, p, beta, flag);
      add_vec(size, p, z, flag);
    } else { 
      copy_vec(size, z, p, flag);
    }

    rsb_SPMV(spMatrix, p, work, flag);

    alpha = rho/ddot_product(size, p, work);
    scale_vec(size, p, alpha, flag);
    add_vec(size, sol, p, flag);

    scale_vec(size, work, alpha, flag);
    sub_vec(size, r, work, flag);
    err = sqrt(ddot_product(size, r, r))/bNorm2;

    if (err <= tol) {
      printf("%d\n", i);
      return NULL;
    }
    rho1 = rho;
  }

  safe_free( r );
  safe_free( z );
  safe_free( p );
  safe_free( work );


   
  return NULL;
} 


void *sp_block_conj_gradd(struct rsb_mtx_t *spMatrix, double *dMat,                       
                          double *solMat, const int size, const int nrhs,                 
                          struct rsb_mtx_t *spMatrixDiag, const int maxIter,                      
                          const double tol) {                                                     
  
  // solves the symmetric positive definite linear system AX=B                                    
  // using the block conjugate gradient method (with preconditioning)                             
                                                                                                  
  double err;                                                                                     
                                                                                                  
  double *alpha, *rho, *beta, *q;                                                         
  alpha = xmalloc(nrhs*nrhs*sizeof(double));                                              
  rho = xmalloc(nrhs*nrhs*sizeof(double));                                                
  q = xmalloc(size*nrhs*sizeof(double));                                                  
  beta = xmalloc(nrhs*nrhs*sizeof(double));                                               
                                                                                                  
  const char flag = DATA_TYPE_DOUBLE;    
                                          
  // get norm of B                        
  // if norm of B == 0, return null vector
  double BNorm2 = matrix_norm(size, nrhs, dMat, 'F', flag);

  if (BNorm2 < 1.0e-30) {
    BNorm2 = 1.0;
    err = 0.0;
    memset(solMat, 0, size*nrhs*sizeof(double));

    return NULL;
  }
                                      
  double *r, *z, *p;          
  r = xmalloc(size*nrhs*sizeof(double));
  z = xmalloc(size*nrhs*sizeof(double));
  p = xmalloc(size*nrhs*sizeof(double));
                                                                                                  
  // why is this here?                                                                            
  copy_vec(size*nrhs, dMat, r, flag);

  double *work, *work1;
  work = xmalloc(nrhs*nrhs*sizeof(double));
  work1 = xmalloc(nrhs*nrhs*sizeof(double));
                                                                                                  
  // calculate residual r                                                                         
  rsb_SPMM_scal_add(spMatrix, solMat, r, nrhs, -1.0, 1.0, flag);                                  
  err = matrix_norm(size, nrhs, r, 'F', flag)/BNorm2;                                             
                                                                                                  
  if(err <= tol)
    return NULL;
  
  for(int i = 0; i < maxIter; i++) {

    // orthogonalize preconditioned residual
    rsb_SPMM(spMatrixDiag, r, z, nrhs, flag);
    if(i > 0) {
      // calculate beta = -(P(i)'*Q(i))^-1(Q(i)'*Z(i+1)) 
      // where (P(i)'*Q(i))^-1 is still in work                                          
      gemm_TN(nrhs, nrhs, size, q, z, work1, flag);
      gemm_NN_scal_add(nrhs, nrhs, nrhs, work, work1, beta, -1.0, 0.0, flag);
      gemm_NN_scal_add(size, nrhs, nrhs, p, beta, z, 1.0, 1.0, flag);        
    }
    
    get_q(size, nrhs, z, flag);
    copy_vec(size*nrhs, z, p, flag);
    
    // calculate alpha as P'*R*(P'Q)^-1
    // where Q = A*P
    gemm_TN(nrhs, nrhs, size, p, r, rho, flag);
    rsb_SPMM(spMatrix, p, q, nrhs, flag);
                                         
    gemm_TN(nrhs, nrhs, size, p, q, work, flag);
    matrix_inverse(nrhs, work, flag);           
    gemm_NN(nrhs, nrhs, nrhs, work, rho, alpha, flag);

    // get new solution X(i+1) = X(i)+P(i)*alpha
    gemm_NN_scal_add(size, nrhs, nrhs, p, alpha, solMat, 1.0, 1.0, flag);
                                                                                                  
    // get new residual R(i+1) = R(i) - Q(i)*alpha
    gemm_NN_scal_add(size, nrhs, nrhs, q, alpha, r, -1.0, 1.0, flag);
    err = matrix_norm(size, nrhs, r, 'F', flag)/BNorm2;
                                                                                                  
    if (err <= tol)
      return NULL; 
                   
  }                

  safe_free( r );
  safe_free( z );
  safe_free( p );
  safe_free( q );
  safe_free( alpha );
  safe_free( beta ); 
  safe_free( rho );
                   
  return NULL;     
} 

void *sp_block_conj_gradz(struct rsb_mtx_t *spMatrix, double complex *zMat,
			  double complex *solMat, const int size, const int nrhs,
			  struct rsb_mtx_t *spMatrixDiag, const int maxIter,
			  const double tol) {

  // soles the symmetric positive definite linear system AX=B
  // using the breakdown-free block conjugate gradient method
  // (with preconditioning)

  double err;

  double complex *alpha, *rho, *beta, *q;
  alpha = xmalloc(nrhs*nrhs*sizeof(double complex));
  rho = xmalloc(nrhs*nrhs*sizeof(double complex));
  q = xmalloc(size*nrhs*sizeof(double complex));
  beta = xmalloc(nrhs*nrhs*sizeof(double complex));

  const char flag = DATA_TYPE_COMPLEX;

  // get norm of B
  // if norm of B == 0, return null vector
  double BNorm2 = matrix_norm(size, nrhs, zMat, 'F', flag);

  if (BNorm2 < 1.0e-30) {
    BNorm2 = 1.0;
    err = 0.0;
    memset(solMat, 0, size*nrhs*sizeof(double complex));

    safe_free( alpha );
    safe_free( rho );
    safe_free( q );
    safe_free( beta );

    return NULL;
  }

  double complex *r, *z, *p;
  r = xmalloc(size*nrhs*sizeof(double complex));
  z = xmalloc(size*nrhs*sizeof(double complex));
  p = xmalloc(size*nrhs*sizeof(double complex));

  copy_vec(size*nrhs, zMat, r, flag);

  double complex *work, *work1;
  work = xmalloc(nrhs*nrhs*sizeof(double complex));
  work1 = xmalloc(nrhs*nrhs*sizeof(double complex));

  // caluclate residual r
  rsb_SPMM_scal_add(spMatrix, solMat, r, nrhs, -1.0, 1.0, flag);
  err = matrix_norm(size, nrhs, r, 'F', flag)/BNorm2;

  if(err <= tol) {
    // could potentially turn this into a error per nrhs
    // instead of a total one
    safe_free( alpha );
    safe_free( beta );
    safe_free( rho );
    safe_free( r );
    safe_free( z );
    safe_free( p );
    safe_free( q );
    safe_free( work );
    safe_free( work1 );
    return NULL;
  }

  for(int i = 0; i < maxIter; i++) {

    // orthogonalize preconditioned residual
    rsb_SPMM(spMatrixDiag, r, z, nrhs, flag);

    if(i > 0) {
      // calculate beta = -(P(i)'*Q(i))^-1*(Q(i)'*Z(i+1))
      // where (P(i)'*Q(i))^-1 is still in work
      gemm_TN(nrhs, nrhs, size, q, z, work1, flag);
      gemm_NN_scal_add(nrhs, nrhs, nrhs, work, work1, beta, -1.0, 0.0, flag);
      gemm_NN_scal_add(size, nrhs, nrhs, p, beta, z, 1.0, 1.0, flag);
    }

    get_q(size, nrhs, z, flag);
    copy_vec(size*nrhs, z, p, flag);

    // calculate alpha as P'*R*(P'Q)^-1
    // where Q = A*P
    gemm_TN(nrhs, nrhs, size, p, r, rho, flag);
    rsb_SPMM(spMatrix, p, q, nrhs, flag);

    gemm_TN(nrhs, nrhs, size, p, q, work, flag);
    matrix_inverse(nrhs, work, flag);
    gemm_NN(nrhs, nrhs, nrhs, work, rho, alpha, flag);

    // get new solution (X+1) = X(i)+P(i)*alpha
    gemm_NN_scal_add(size, nrhs, nrhs, p, alpha, solMat, 1.0, 1.0, flag);

    // get new residual R(i+1) = R(i) - Q(i)*alpha
    gemm_NN_scal_add(size, nrhs, nrhs, q, alpha, r, -1.0, 1.0, flag);
    err = matrix_norm(size, nrhs, r, 'F', flag)/BNorm2;

    if(err <= tol) {
      safe_free( alpha );
      safe_free( beta );
      safe_free( rho );
      safe_free( r );
      safe_free( z );
      safe_free( p );
      safe_free( q );
      safe_free( work );
      safe_free( work1 );
      return NULL;
    }
  }

  safe_free( alpha );
  safe_free( beta );
  safe_free( rho );
  safe_free( r );
  safe_free( z );
  safe_free( p );
  safe_free( q );
  safe_free( work );
  safe_free( work1 );

  return NULL;
}
