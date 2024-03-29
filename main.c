#include "./src/splrep.h"
#include "./src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  const int type = 2;
  const double shift = 0.0;
  const int rsbTune = 0;

  const int dim = 1;

  int N[dim];
  N[0] = 256;
  N[1] = 128;

  double ln[dim];
  ln[0] = 20;
  ln[1] = ln[0];

  const int nev = 20;

  const double mu = 1.4134287173e+01;//5.7597536547e+00;//1.7888646026e+01;
  double param[5] = { 0.0 };
  param[0] = 1.0;
  int trapChoice = 1;
  int intChoice = 1;
  double intParam[5] = { 0.0 };
  intParam[0] = 100.0;

  const double tol = 1.0e-6;
  const int maxIter = 100;

  struct sp_lrep_t *LREP;
  printf("setting up LREP.\n");
  LREP = splrep_setup("filepsi_1D_0.dat", N, ln, nev, dim, mu, param,
		      type, shift, maxIter, tol, rsbTune, DATA_TYPE_DOUBLE);

  printf("setup LREP.\n");

  sp_lopb4dcg(LREP);

  
  /*  //  struct rsb_mtx_t *mtxAp = NULL; // matrix structure pointer
  rsb_coo_idx_t nrA = LREP->size; // number of rows
  rsb_coo_idx_t ncA = LREP->size; // number of cols

  // spmv specific variables
  const RSB_DEFAULT_TYPE alpha = 1;
  const RSB_DEFAULT_TYPE beta = 0;
  rsb_flags_t order = RSB_FLAG_WANT_COLUMN_MAJOR_ORDER;
  const rsb_coo_idx_t nrhs = nev; // number of right hand sides
  rsb_trans_t transA = RSB_TRANSPOSITION_N; // transposition
  rsb_nnz_idx_t ldB = nrA;
  rsb_nnz_idx_t ldC = ncA;

  // misc variables
  rsb_err_t errval = RSB_ERR_NO_ERROR;
  rsb_time_t dt, odt;
  rsb_int_t t, tt = 100; // will repeat spmv tt times
  char ib[200];
  const char *is = "RSB_MIF_MATRIX_INFO__TO__CHAR_P";

  // input autotuning variables
  rsb_int_t oitmax = 100; // autotune iterations
  rsb_time_t tmax = 1.0; // time per autotune operation

  // input/output autotuning variables
  rsb_int_t tn = 0; // threads number

  // output autotuning variables
  rsb_real_t sf = 0.0; // speedup factor obtained from auto tuning
  rsb_int_t wvat = 1; // want verbose autotuning

  errval = rsb_lib_set_opt(RSB_IO_WANT_VERBOSE_TUNING, &wvat);
  //  mtxAp = LREP->K;

  dt = -rsb_time();
  for(t = 0; t < tt; ++t)
    rsb_SPMM(LREP->K, LREP->X, LREP->Y, nev, 'D');
  dt += rsb_time();
  odt = dt;
  printf("Before auto-tuning, %d multiplications took %lfs,\n", tt, dt);
  printf("Threads autotuning (may take more than %lfs) ... \n",
	 oitmax*tmax);

  dt = -rsb_time();
  errval = rsb_tune_spmm(NULL, &sf, &tn, oitmax, tmax, transA,
			 &alpha, LREP->K, nrhs, order, LREP->X, ldB, &beta,
			 LREP->Y, ldC);
  dt += rsb_time();
  if(errval != RSB_ERR_NO_ERROR)
    goto err;

  if(tn == 0)
    printf("After %lfs, autotuning routine did not find a better"
	   "threads count configuration.\n", dt);
  else
    printf("After %lfs, autotuning routine declared speedup of %lg x,"
	   " when using threads count of %d.\n", dt, sf, tn);

  errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn);
  if(errval != RSB_ERR_NO_ERROR)
    goto err;

  rsb_mtx_get_info_str(LREP->K, is, ib, sizeof(ib));
  printf("%s\n", ib);

  dt = -rsb_time();
  for(t = 0; t<tt; ++t)
    rsb_SPMM(LREP->K, LREP->X, LREP->Y, nev, 'D');
  dt += rsb_time();
  printf("\n After threads auto-tuning, %d multiplication took %lfs"
	 " -- effective speedup of %lg x\n", tt, dt, odt/dt);
  
  odt = dt;
  
  tn = 0; // this will restore default threads count
  errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn);
  if(errval != RSB_ERR_NO_ERROR)
    goto err;
  errval = rsb_lib_get_opt(RSB_IO_WANT_EXECUTING_THREADS, &tn);
  if(errval != RSB_ERR_NO_ERROR)
    goto err;

  printf("Matrix autotuning (may take more than %lfs; using %d"
	 "threads ) ...\n", oitmax*tmax, tn);
  
  // A negative tn will request also threads autotuning:
  // tn = -tn;

  dt = rsb_time();
  errval = rsb_tune_spmm(&(LREP->K), &sf, NULL, oitmax, tmax, transA,
			 &alpha, NULL, nrhs, order, LREP->X,
			 ldB, &beta, LREP->Y, ldC);
  dt += rsb_time();

  if(errval != RSB_ERR_NO_ERROR)
    goto err;

  if(tn == 0)
    printf("After %lfs, autotuning routine did not find a better"
	   " structure.\n", dt);
  else
    printf("After %lfs, autotuning routine declared speedup of %lg x,"
	   " when using threads count of %d.\n", dt,sf,tn);

  rsb_mtx_get_info_str(LREP->K, is, ib, sizeof(ib));
  printf("%s\n", ib);

  dt = -rsb_time();
  for(t=0; t<tt; ++t)
    rsb_SPMM(LREP->K, LREP->X, LREP->Y, nev, 'D');
  dt += rsb_time();
  printf("\n After threads auto-tuning, %d multiplications took %lfs"
	 " -- further speed up of %lg x\n", tt, dt, odt/dt);

  //  rsb_mtx_free(mtxAp);
  */
  splrep_free( LREP );
  return 0;

  //err:
  //  rsb_perror(NULL,errval);
  //  printf("Program terminating with error.\n");
  //  return -1;
  
  /*
  // Trying to test all my complex blas and lapack functions
  double complex a[9];
  a[0] = 1.0 + 9.0*I;
  a[1] = 2.0 + 8.0*I;
  a[2] = 3.0 + 7.0*I;

  a[3] = 4.0 + 6.0*I;
  a[4] = 5.0 + 5.0*I;
  a[5] = 6.0 + 4.0*I;

  a[6] = 7.0 + 3.0*I;
  a[7] = 8.0 + 2.0*I;
  a[8] = 9.0 + 1.0*I;

  double complex b[9];
  b[0] = 1.0 - I;
  b[1] = 0.0;
  b[2] = 0.0;
  
  b[3] = 0.0;
  b[4] = 1.0 - 2*I;
  b[5] = 0.0;
  
  b[6] = 0.0;
  b[7] = 0.0;
  b[8] = 1.0 - 3*I;

  double complex out[9] = { 0.0 + I*0.0 };

  double dtest[9];
  dtest[0] = 1.0;
  dtest[1] = 1.0;
  dtest[2] = 1.0;
  dtest[3] = 1.0;
  dtest[4] = 1.0;
  dtest[5] = 1.0;
  dtest[6] = 1.0;
  dtest[7] = -8.0;
  dtest[8] = 8.0;

  double ret = 0.0;
  ret = dasum(9, a, DATA_TYPE_COMPLEX);
  printf("%.5e\n", ret);
  */
  return 0;
}
