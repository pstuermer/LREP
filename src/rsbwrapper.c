
#include "rsbwrapper.h"
#include <rsb.h>

void rsb_init(const int rsbTune) {
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  errVal = rsb_lib_init(RSB_NULL_INIT_OPTIONS);
  if(errVal != RSB_ERR_NO_ERROR)
    printf("Error initializing rsb_lib: %d", errVal);

  if(rsbTune == 1) {
    errVal = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &rsbTune);
    if(errVal != RSB_ERR_NO_ERROR)
      printf("Error setting verbose tuning: %d", errVal);
  }
}

struct rsb_mtx_t *rsb_mtx_from_coo(coo_t *matrix) {
  struct rsb_mtx_t* mtxAp = NULL;

  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t dtypeCode = RSB_NUMERICAL_TYPE_DEFAULT;
  rsb_type_t ztypeCode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX;

  if(matrix->flag == 'D') {
    mtxAp = rsb_mtx_alloc_from_coo_const(matrix->dval, matrix->row, matrix->col,
					 matrix->nnz, dtypeCode, matrix->rows,
					 matrix->cols, brA, bcA,
					 RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS, &errVal);
  } else if(matrix->flag == 'Z') {
    mtxAp = rsb_mtx_alloc_from_coo_const(matrix->zval, matrix->row, matrix->col,
					 matrix->nnz, ztypeCode, matrix->rows,
					 matrix->cols, brA, bcA,
					 RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS, &errVal);
  }

  if((!mtxAp) || (errVal != RSB_ERR_NO_ERROR)) {
    printf("Error allocating RSB matrix from coo: %d", errVal);
    rsb_perror(NULL, errVal);
  }

  return mtxAp;
}

struct rsb_mtx_t *rsb_mtx_from_coo_sym(coo_t *matrix) {
  struct rsb_mtx_t *mtxAp = NULL;

  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t dtypeCode = RSB_NUMERICAL_TYPE_DEFAULT;
  rsb_type_t ztypeCode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX;

  if(matrix->flag == 'D')
    mtxAp = rsb_mtx_alloc_from_coo_const(matrix->dval, matrix->row, matrix->col,
					 matrix->nnz, dtypeCode, matrix->rows,
					 matrix->cols, brA, bcA,
					 RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
					 | RSB_FLAG_DUPLICATES_SUM
					 | RSB_FLAG_DISCARD_ZEROS
					 | RSB_FLAG_SYMMETRIC,
					 &errVal);
  if(matrix->flag == 'Z')
    mtxAp = rsb_mtx_alloc_from_coo_const(matrix->dval, matrix->row, matrix->col,
					 matrix->nnz, ztypeCode, matrix->rows,
					 matrix->cols, brA, bcA,
					 RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
					 | RSB_FLAG_DUPLICATES_SUM
					 | RSB_FLAG_DISCARD_ZEROS
					 | RSB_FLAG_SYMMETRIC,
					 &errVal);
  if((!mtxAp) || (errVal != RSB_ERR_NO_ERROR)) {
    printf("Error allocating symmetric RSB matrix from coo: %d", errVal);
    rsb_perror(NULL, errVal);
  }

  return mtxAp;
}
  

void rsb_SPMM(struct rsb_mtx_t *spMatrix, const void *dMatrix,
	      void *res, int nev, const char flag) {

  const RSB_DEFAULT_TYPE zero = 0;
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zzero = 0;
  const double complex zone = 1;
  
  rsb_err_t errVal = RSB_ERR_NO_ERROR;

  if(flag == 'D') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zero,res,nev);

  } else if(flag == 'Z') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zzero,res,nev);
  }

  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPMM: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}
  
void rsb_SPMM_sub(struct rsb_mtx_t *spMatrix, const void *dMatrix,
		  void *res, int nev, const char flag) {

  const RSB_DEFAULT_TYPE minus = -1;  
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zminus = -1;
  const double complex zone = 1;
 
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  
  if(flag == 'D') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &minus,res,nev);
  } else if(flag == 'Z') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zminus,res,nev);
  }
  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPMM_sub: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}

void rsb_SPMM_scal_add(struct rsb_mtx_t *spMatrix, const void *dMatrix,
		       void *res, int nev, const double alpha,
		       const double beta, const char flag) {

  const RSB_DEFAULT_TYPE dalpha = alpha;
  const RSB_DEFAULT_TYPE dbeta = beta;
  const double complex zalpha = alpha;
  const double complex zbeta = beta;  
   
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
   
  if(flag == 'D') {    
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&dalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &dbeta,res,nev);
  } else if(flag == 'Z') {  
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &zbeta,res,nev);
  }
  if(errVal != RSB_ERR_NO_ERROR) {    
    printf("Error performing a SPMM_scal_add: %d\n", errVal);  
    rsb_perror(NULL,errVal);
  }
}
   
void rsb_SPMV(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag) {   
   
  const RSB_DEFAULT_TYPE zero = 0;    
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zzero = 0;
  const double complex zone = 1;
  
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
 
  if(flag == 'D') {
    errVal = rsb_spmv(RSB_TRANSPOSITION_N,&one,spMatrix,vec,1,
		      &zero,res,1);
  } else if (flag == 'Z') {
    errVal = rsb_spmv(RSB_TRANSPOSITION_N,&zone,spMatrix,vec,1,
		      &zzero,res,1);
  }

  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPMV: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}

void rsb_SPVM(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag) {
   
  const RSB_DEFAULT_TYPE zero = 0;
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zzero = 0;
  const double complex zone = 1;
  
  rsb_err_t errVal = RSB_ERR_NO_ERROR;

  if(flag == 'D') {
    errVal = rsb_spmv(RSB_TRANSPOSITION_T,&one,spMatrix,vec,1,
		      &zero,res,1);
  } else if (flag == 'Z') {
    errVal = rsb_spmv(RSB_TRANSPOSITION_C,&zone,spMatrix,vec,1,
		      &zzero,res,1);
  }

  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPVM: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}  
   
void rsb_SPSM(struct rsb_mtx_t *spMatrix, const void *dMatrix, void *res,
	      int nev, const char flag) {    
  const RSB_DEFAULT_TYPE zero = 0;    
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zzero = 0;
  const double complex zone = 1;
  
  rsb_err_t errVal = RSB_ERR_NO_ERROR;

  if(flag == 'D') {
    errVal = rsb_spsm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,&zero,dMatrix,
		      nev,res,nev);
  } else if(flag == 'Z') {
    errVal = rsb_spsm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,&zzero,dMatrix,
		      nev,res,nev);
  }

  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPSM: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}  
   
void rsb_SPSV(struct rsb_mtx_t *spMatrix, const void *vec, void *res,
	      const char flag) {   
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zone = 1;
  
  rsb_err_t errVal = RSB_ERR_NO_ERROR;

  if(flag == 'D') {
    rsb_spsv(RSB_TRANSPOSITION_N,&one,spMatrix,vec,1,res,1);
  } else if (flag == 'Z') {
    rsb_spsv(RSB_TRANSPOSITION_N,&zone,spMatrix,vec,1,res,1);
  }

  if(errVal != RSB_ERR_NO_ERROR) {
    printf("Error performing a SPSV: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}
   
void rsb_get_prec(struct rsb_mtx_t *spMatrix, void *opdp[2]) {
  const void *ipdp[2]; 
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
   
  if(( errVal = rsb_mtx_get_prec(opdp, spMatrix, RSB_PRECF_ILU0, ipdp)) !=   
     RSB_ERR_NO_ERROR) {    
    printf("Error calculating a preconditioner: %d\n", errVal);
    rsb_perror(NULL,errVal);
  }
}

void rsb_tune_SPMM(struct rsb_mtx_t *spMatrix, const double *dMatrix, int nev,    
		   int size, int tn, const char flag) {
  int tt = 100;   
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const RSB_DEFAULT_TYPE zero = 0;    
  const RSB_DEFAULT_TYPE one = 1;
  rsb_trans_t transA = RSB_TRANSPOSITION_N;
  rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER; 
  rsb_time_t dt;  
  rsb_time_t dtpre;    
  const rsb_int_t oitmax = 100; // number of autotune iterations
  const rsb_time_t tmax = 0.1; // time per autotune iteration in s  
  //  struct rsb_mtx_t* mtxOp = NULL; // new matrix pointer for better matrix structure
  rsb_real_t sf = 0.0; // ideal speed-up factor   
  rsb_int_t tn1 = 0;   
  char ib[200];   
   
  double *res;    
  res = xmalloc((size*nev+1)*sizeof(double));
   
  // Do one matrix multiplication for actual time value   
  dtpre = -rsb_time(); 
  for(int t = 0; t < tt; t++)    
    rsb_SPMM(spMatrix, dMatrix, res, nev, flag);    
  dtpre += rsb_time(); 
   
  printf("Before auto-tuning, %d multiplications sparse*dense took %lfs.\n",tt,dtpre);  
  printf("Threads autotuning (may take more than %lfs)...\n", oitmax*tmax);   
  
  dt = -rsb_time();
  errVal = rsb_tune_spmm(NULL, &sf, &tn1, oitmax, tmax, transA, &one,    
			 spMatrix, nev, order, dMatrix, size, &zero,    
			 res, size);   
  dt += rsb_time();    
  if(errVal != RSB_ERR_NO_ERROR) 
    goto err;
   
  if(tn == 0)
    printf("After %lfs, autotuning routine did not find a better"   
	   " threads count configuration.\n",dt);    
  else  
    printf("After %lfs, autotuning routine declared speedup of %lg x,"   
	   " when using threads count of %d.\n",dt,sf,tn);
    
  errVal = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS,&tn); 
  if(errVal != RSB_ERR_NO_ERROR) 
    goto err;
   
  dt = -rsb_time();    
  for(int t = 0; t < tt; t++)    
    rsb_SPMM(spMatrix, dMatrix, res, nev, flag);    
  dt += rsb_time();    
  printf("After threads auto-tuning, %d multiplications took %lfs"  
	 "  --  effective speedup of %lg x\n",tt,dt,dtpre/dt); 
  dtpre = dt;
   
  // restore default thread counts    
  errVal = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn1);    
  if(errVal != RSB_ERR_NO_ERROR) 
    goto err;
  errVal = rsb_lib_get_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn1);    
  if(errVal != RSB_ERR_NO_ERROR) 
    goto err;
   
   
  printf("Matrix autotuning (may take more than %lfs; using %d"
	 " threads )...\n", oitmax*tmax, tn1);  
  //  mtxOp = spMatrix;
  dt = -rsb_time();

  errVal = rsb_tune_spmm(&spMatrix,&sf,&tn,oitmax,tmax,transA,
			 &one,NULL,nev,order,NULL,size,&zero,
			 NULL,size);
  
  dt += rsb_time();
  if(errVal != RSB_ERR_NO_ERROR) 
    goto err;
   
  /*  if( mtxOp == spMatrix )  
      {   
      printf("After %lfs, global autotuning found old matrix optimal,"   
      " with declared speedup %lg x when using %d threads\n",dt,sf,tn);
      }   
      else  
      {   
      printf("After %lfs, global autotuning declared speedup of %lg x,"  
      " when using threads count of %d and a new matrix:\n",dt,sf,tn); 
      }   
  */
  errVal = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn);
   
  dt = -rsb_time();    
  for(int t = 0; t < tt; t++)    
    rsb_SPMM(spMatrix, dMatrix, res, nev, flag);    
  dt += rsb_time();    
   
  printf("After global auto-tuning, %d multiplications took %lfs"   
	 "  --  further speedup of %lg x\n",tt,dt,dtpre/dt);   
   
  // restore default thread counts    
  //errVal = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn1);    
  //if(errVal != RSB_ERR_NO_ERROR) 
  // goto err;
   
  //  rsb_mtx_free( mtxOp );    
  safe_free( res ); 
   
 err:   
  rsb_perror(NULL, errVal); 
  rsb_strerror_r(errVal, ib, 200);    
  printf("Error tuning a SPMM: %d\n", errVal);  
  printf("rsb_sterror_r: %s\n", ib);  
  printf("Program terminating with error.\n");  
  exit(1);   
}

/*
  struct rsb_mtx_t *rsb_mtx_from_zcoo(zcoo_t *matrix) {
  struct rsb_mtx_t *mtxAp = NULL;

  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t typeCode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX;

  mtxAp = rsb_mtx_alloc_from_coo_const(matrix->val, matrix->row, matrix->col,
  matrix->nnz, typeCode, matrix->rows,
  matrix->cols, brA, bcA, RSB_FLAG_NOFLAGS,
  &errVal);

  if((!mtxAp) || (errVal != RSB_ERR_NO_ERROR)) {
  printf("Error allocating RSB matrix from coo: %d\n", errVal);
  rsb_perror(NULL, errVal);
  }
  return mtxAp;
  }
*/
