#include "splrep.h"
#include "precond.h"
#include <rsb.h>

struct sp_lrep_t *splrep_malloc(const int *N, const int nev, const int sizeSub,
				const int dim, coo_t *K, coo_t *M, 
				const int type, const double shift,
				const int maxIter, const double tol,
				const char flag) {
  
  struct sp_lrep_t *LREP = xmalloc(sizeof(struct sp_lrep_t));
  
  LREP -> dim = dim;
  LREP -> nev = nev;
  LREP -> sizeSub = sizeSub;
  LREP -> iter = 0;
  LREP -> maxIter = maxIter;
  LREP -> tol = tol;
  LREP -> flag = flag;
  
  int size = 1;
  LREP -> N = xmalloc(dim*sizeof(int));
  for(int i = 0; i < dim; i++) {
    LREP -> N[i] = N[i];
    size *= N[i];
  }
  LREP -> size = size;
  
  if (LREP->flag == 'Z') {
  
    LREP -> X = xmalloc(size*sizeSub*sizeof(double complex));
    LREP -> Y = xmalloc(size*sizeSub*sizeof(double complex));
    init_eig_vec(size, sizeSub, LREP->X, -1.0, 1.0, flag);
    init_eig_vec(size, sizeSub, LREP->Y, -1.0, 1.0, flag);
    LREP -> X1 = xmalloc(size*sizeSub*sizeof(double complex));
    LREP -> Y1 = xmalloc(size*sizeSub*sizeof(double complex));
    
    LREP -> eigVal = calloc(sizeSub*sizeSub, sizeof(double complex));
    
    LREP -> P = xmalloc(size*sizeSub*sizeof(double complex));
    LREP -> Q = xmalloc(size*sizeSub*sizeof(double complex));
    
    LREP -> U = xmalloc(3*size*sizeSub*sizeof(double complex));
    LREP -> V = xmalloc(3*size*sizeSub*sizeof(double complex));
    
    LREP -> W = xmalloc(3*3*sizeSub*sizeSub*sizeof(double complex));
    LREP -> Hsr = xmalloc(6*6*sizeSub*sizeSub*sizeof(double complex));
    
    LREP -> eVecsr = xmalloc(6*6*sizeSub*sizeSub*sizeof(double complex));
    LREP -> eVecSort = xmalloc(6*sizeSub*sizeSub*sizeof(double complex));
    LREP -> Xsr = xmalloc(3*sizeSub*sizeSub*sizeof(double complex));
    LREP -> Ysr = xmalloc(3*sizeSub*sizeSub*sizeof(double complex));
    LREP -> eValsr = xmalloc(6*sizeSub*sizeof(double complex));
    LREP -> eValSort = xmalloc(6*sizeSub*sizeof(double complex));
  }
  LREP -> resNorm = xmalloc(sizeSub*sizeof(double));

  LREP -> K = rsb_mtx_from_coo_sym(K);
  LREP -> M = rsb_mtx_from_coo_sym(M);
  
  LREP -> spPrecond = cond_malloc(type, shift);
  
  double KNorm;
  double MNorm;
  rsb_mtx_get_norm(LREP->K, &KNorm, RSB_EXTF_NORM_ONE);
  rsb_mtx_get_norm(LREP->M, &MNorm, RSB_EXTF_NORM_ONE);
  
  LREP -> oneNorm = MAX(KNorm, MNorm);
  
  return LREP;
}

struct sp_lrep_t *rsb_splrep_malloc(const int *N, const int nev, const int sizeSub,
				    const int dim, struct rsb_mtx_t *K, struct rsb_mtx_t *M, 
				    const int cond_type, const double shift,
				    const int maxIter, const double tol,
				    const char flag) {
  
  struct sp_lrep_t *LREP = xmalloc(sizeof(struct sp_lrep_t));

  LREP -> dim = dim;
  LREP -> nev = nev;
  LREP -> sizeSub = sizeSub;
  LREP -> iter = 0;
  LREP -> maxIter = maxIter;
  LREP -> tol = tol;
  LREP -> flag = flag;
  
  int size = 1;
  LREP -> N = xmalloc(dim*sizeof(int));
  for(int i = 0; i < dim; i++) {
    LREP -> N[i] = N[i];
    size *= N[i];
  }
  LREP -> size = size;

  
  LREP -> X = xmalloc(size*sizeSub*sizeof(double));
  LREP -> Y = xmalloc(size*sizeSub*sizeof(double));
  init_eig_vec(size, sizeSub, LREP->X, -1.0, 1.0, flag);
  init_eig_vec(size, sizeSub, LREP->Y, -1.0, 1.0, flag);
  LREP -> X1 = xmalloc(size*sizeSub*sizeof(double));
  LREP -> Y1 = xmalloc(size*sizeSub*sizeof(double));
    
  LREP -> eigVal = calloc(sizeSub*sizeSub, sizeof(double));
    
  LREP -> P = xmalloc(size*sizeSub*sizeof(double));
  LREP -> Q = xmalloc(size*sizeSub*sizeof(double));
    
  LREP -> U = xmalloc(3*size*sizeSub*sizeof(double));
  LREP -> V = xmalloc(3*size*sizeSub*sizeof(double));
    
  LREP -> W = xmalloc(3*3*sizeSub*sizeSub*sizeof(double));
  LREP -> Hsr = xmalloc(6*6*sizeSub*sizeSub*sizeof(double));
    
  LREP -> eVecsr = xmalloc(6*6*sizeSub*sizeSub*sizeof(double));
  LREP -> eVecSort = xmalloc(6*sizeSub*sizeSub*sizeof(double));
  LREP -> Xsr = xmalloc(3*sizeSub*sizeSub*sizeof(double));
  LREP -> Ysr = xmalloc(3*sizeSub*sizeSub*sizeof(double));
  LREP -> eValsr = xmalloc(6*sizeSub*sizeof(double));
  LREP -> eValSort = xmalloc(6*sizeSub*sizeof(double));
  LREP -> resNorm = xmalloc(sizeSub*sizeof(double));

  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  rsb_type_t rsb_type = RSB_NUMERICAL_TYPE_SAME_TYPE;
  rsb_flags_t flags = RSB_FLAG_IDENTICAL_FLAGS;
  rsb_trans_t trans = RSB_TRANSPOSITION_N;
  errVal = rsb_mtx_clone(&(LREP->K), rsb_type, trans, NULL, K, flags);
  if( RSB_ERR_NO_ERROR != errVal )
    printf("Problem cloning K matrix when allocating LREP.\n");
  errVal = rsb_mtx_clone(&(LREP->M), rsb_type, trans, NULL, M, flags);
  if( RSB_ERR_NO_ERROR != errVal )
    printf("Problem cloning M matrix when allocating LREP.\n");
  
  LREP -> spPrecond = cond_malloc(cond_type, shift);
  
  double KNorm;
  double MNorm;
  rsb_mtx_get_norm(LREP->K, &KNorm, RSB_EXTF_NORM_ONE);
  rsb_mtx_get_norm(LREP->M, &MNorm, RSB_EXTF_NORM_ONE);
  
  LREP -> oneNorm = MAX(KNorm, MNorm);
  
  return LREP;
}

void splrep_free(sp_lrep_t *LREP){ 
  assert(LREP != NULL); 

  //  printf("1\n");
  
  rsb_mtx_free( LREP -> K );
  rsb_mtx_free( LREP -> M );
  safe_free( LREP -> N );
  //  printf("2\n");
  
  safe_free( LREP -> X );
  safe_free( LREP -> Y );
  safe_free( LREP -> X1 );
  safe_free( LREP -> Y1 );
  safe_free( LREP -> P );
  safe_free( LREP -> Q );
  safe_free( LREP -> U );
  safe_free( LREP -> V );
  safe_free( LREP -> W );
  safe_free( LREP -> Hsr );
  //  printf("3\n");
  
  safe_free( LREP -> eVecsr );
  //  printf("1\n");
  safe_free( LREP -> eVecSort ); 
  //  printf("2\n");
  safe_free( LREP -> Xsr );
  //  printf("3\n");
  safe_free( LREP -> Ysr );
  //  printf("4\n");
  safe_free( LREP -> eValsr );
  //  printf("5\n");
  safe_free( LREP -> eValSort ); 
  //  printf("6\n");
  safe_free( LREP -> resNorm );
  //  printf("4\n");
  
  safe_free( LREP -> eigVal );
  cond_free( LREP -> spPrecond );
  //  printf("5\n");
  
  safe_free( LREP ); 
  //  printf("6\n");
  rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
}

struct sp_lrep_t *splrep_setup(char fileName[], const int *N, const double *ln,
			       const int nev, const int subspaceSize,
			       const int dim, const double mu,
			       const int trapChoice, const double *param,
			       const int intChoice, const double *intParam,
			       const int type, const double shift,
			       const int maxIter, const double tol,
			       const int rsbTune, const char flag) {  
  // load wf
  double complex *wf;
  double *apsi; 
  int factor = 1, size = 1; 
  int NFac[dim];

  for(int i = 0; i < dim; i++) {
    size *= N[i];
    NFac[i] = N[i]/factor; 
  }

  wf = xmalloc(size/factor*sizeof(double complex));
  
  load_wf(wf, size, factor, fileName); 
 
  size /= factor;
  apsi = get_apsi(wf, size); 

  // set-up grid
  grid_t *grid; 
  grid = grid_malloc(dim, NFac, ln); 
  grid_setup(grid); 
  
  // Setup differentiatior in COO format
  coo_t *diff;
  diff = setup_diff2(ln, NFac, dim, flag);

  // so I got until here with generalizing it for 1D and 3D 
  // the question now is, how would I continue for the interaction  
  // and trap parameter 
  // Setup trap and interaction term in COO format  
  coo_t *cooTrap, *cooUint, *cooB1;
  cooTrap = coo_malloc(size, size, size, flag);
  cooUint = coo_malloc(size, size, size, flag);
  cooB1 = coo_malloc(size, size, size, flag);

  setup_trap(cooTrap, NFac, param, grid, dim, trapChoice);
  setup_Uint(cooUint, apsi, NULL, NFac, dim, mu, intParam, intChoice,
	     SYS_ONE_COMPONENT, 0);
  setup_Bint(cooB1, apsi, NULL, NFac, dim, intParam, intChoice,
	     SYS_ONE_COMPONENT, 0);
  
  coo_t *K, *M;
  K = coo_malloc(diff->nnz, size, size, flag);
  M = coo_malloc(diff->nnz, size, size, flag);
  
  setup_KM1C(K, diff, cooTrap, cooUint, cooB1, 'K', size, flag);
  setup_KM1C(M, diff, cooTrap, cooUint, cooB1, 'M', size, flag);
  printf("M done.\n");
  
  // turn K and M from coo into rsb 
  rsb_init(rsbTune);
  struct sp_lrep_t *LREP;
  LREP = splrep_malloc(NFac, nev, subspaceSize, dim, K, M,
		       type, shift, maxIter, tol, flag); 
  
  sp_setup_precond(LREP); 
  printf("setup preconditioner.\n");
  /*  
  // put tuning and optimization into here  
  if(rsbTune == 1) {   
    rsb_int_t tnK = 0;
    rsb_int_t tnM = 0;
    rsb_tune_SPMM(LREP->K, LREP->X, subspaceSize, size, tnK, flag);
    rsb_tune_SPMM(LREP->M, LREP->X, subspaceSize, size, tnM, flag);
    if(LREP->spPrecond->MDiag != NULL) {   
      rsb_int_t tnDiag = 0; 
      rsb_tune_SPMM(LREP->spPrecond->MDiag, LREP->X, subspaceSize,
		    size, tnDiag, flag);
    }   
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tnK);
    if(errval != RSB_ERR_NO_ERROR) {
      rsb_perror(NULL, errval); 
      printf("Program terminating with error.\n");
    }   
  } 
  */
  printf("setup LREP.\n");
  coo_free( M );
  coo_free( K );
  coo_free( cooB1 );
  coo_free( cooUint ); 
  coo_free( cooTrap ); 
  coo_free( diff );

  grid_free( grid );

  safe_free( apsi );
  safe_free( wf );

  return LREP;
}

struct sp_lrep_t *rsb_splrep_setup(char fileName[], const int *N, const double *ln,
				   const int nev, const int subspaceSize,
				   const int dim, const double mu,
				   const double *intParam,
				   double (*V_trap_func)(double*),
				   const int type, const double shift,
				   const int maxIter, const double tol,
				   const char flag) {
  //				   const int rsbTune, const char flag) {  
  // load wf
  double complex *wf = NULL;
  double *apsi = NULL; 
  int factor = 1, size = 1; 
  int NFac[dim];

  for(int i = 0; i < dim; i++) {
    size *= N[i];
    NFac[i] = N[i]/factor; 
  }

#pragma omp parallel
  wf = xmalloc(size/factor*sizeof(double complex));
  
  load_wf(wf, size, factor, fileName); 
 
  size /= factor;
  apsi = get_apsi(wf, size);  
  
  // Setup differentiatior in rsb format
  rsb_lib_init(NULL);
  struct rsb_mtx_t *diff;
  diff = rsb_setup_diff2(ln, NFac, dim);

  // set-up grid
  grid_t *grid;
  grid = grid_malloc(dim, NFac, ln);
  grid_setup(grid);

  // Setup trap and interaction term in rsb format  
  struct rsb_mtx_t *trap = NULL, *uint = NULL, *b1 = NULL;
  double mu_eps = mu + 1.0e-13;
  trap = rsb_setup_trap(V_trap_func, grid);
  uint = rsb_setup_int(apsi, mu_eps, N, dim, intParam);
  b1 = rsb_setup_int(apsi, 0.0, N, dim, intParam);

  // Setup K and M matrices in rsb format
  trap = rsb_dadd(trap, uint);

  struct rsb_mtx_t *K = NULL, *M = NULL;
  K = rsb_dsub(trap, b1);
  K = rsb_dadd(diff, K);

  M = rsb_dadd(trap, b1);
  M = rsb_dadd(diff, M);
    
  // setup final LREP 
  struct sp_lrep_t *LREP;
  LREP = rsb_splrep_malloc(NFac, nev, subspaceSize, dim, K, M,
			   type, shift, maxIter, tol, flag); 
  
  sp_setup_precond(LREP); 
  printf("setup preconditioner.\n");
  /*  
  // put tuning and optimization into here  
  if(rsbTune == 1) {   
    rsb_int_t tnK = 0;
    rsb_int_t tnM = 0;
    rsb_tune_SPMM(LREP->K, LREP->X, subspaceSize, size, tnK, flag);
    rsb_tune_SPMM(LREP->M, LREP->X, subspaceSize, size, tnM, flag);
    if(LREP->spPrecond->MDiag != NULL) {   
      rsb_int_t tnDiag = 0; 
      rsb_tune_SPMM(LREP->spPrecond->MDiag, LREP->X, subspaceSize,
		    size, tnDiag, flag);
    }   
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tnK);
    if(errval != RSB_ERR_NO_ERROR) {
      rsb_perror(NULL, errval); 
      printf("Program terminating with error.\n");
    }   
  } 
  */
  printf("setup LREP.\n");
  rsb_mtx_free( M );
  rsb_mtx_free( K );
  rsb_mtx_free( b1 );
  rsb_mtx_free( uint ); 
  rsb_mtx_free( trap );
  // need to figure out a way to do this if I want to save
  // the differential operator when I do a sweep
  rsb_mtx_free( diff );

  grid_free( grid );

  safe_free( apsi );
  safe_free( wf );

  return LREP;
}

struct sp_lrep_t *splrep2C_setup(char fileName1[], char fileName2[],
				 const int *N, const double *ln,
				 const int nev, const int subspaceSize,
				 const int dim,const double mu[2],
				 const int trapChoice, const double *param,
				 const int intChoice, const double *intParam,
				 const int type, const double shift,
				 const int maxIter, const double tol,
				 const int rsbTune, const char flag) {

  // load wfs
  double complex *wf1, *wf2;
  double *apsi1, *apsi2;

  int size = 1;
  int factor = 1;
  int N2[dim];
  int NFac[dim];
  printf("start.\n");
  for(int i = 0; i < dim; i++) {
    size *= N[i];
    N2[i] = 2*N[i]/factor;
    NFac[i] = N[i]/factor;
  }
  wf1 = xmalloc(size/factor*sizeof(double complex));
  wf2 = xmalloc(size/factor*sizeof(double complex));

  if (dim == 1) {
    load_wf_2C(wf1, wf2, size, factor, fileName1);
  } else {
    load_wf(wf1, size, factor, fileName1);
    load_wf(wf2, size, factor, fileName2);
  }

  size /= factor;
  //  printf("load wf.\n");  
  apsi1 = get_apsi(wf1, size);
  apsi2 = get_apsi(wf2, size);
  //  printf("got apsi.\n");

  // set up grid
  grid_t *grid;
  grid = grid_malloc(dim, NFac, ln);
  grid_setup(grid);

  //  printf("setup grid.\n");

  // set up differentiator in COO format
  coo_t *diff;
  diff = setup_diff2(ln, NFac, dim, flag);

  //  printf("setup diff.\n");

  // set up trap and interaction term for each component in
  // COO format
  coo_t *cooTrap, *cooC, *cooUint1, *cooUint2, *cooB1, *cooB2;
  cooTrap = coo_malloc(size, size, size, flag);
  cooC = coo_malloc(size, size, size, flag);
  cooUint1 = coo_malloc(size, size, size, flag);
  cooUint2 = coo_malloc(size, size, size, flag);
  cooB1 = coo_malloc(size, size, size, flag);
  cooB2 = coo_malloc(size, size, size, flag);

  setup_trap(cooTrap, NFac, param, grid, dim, trapChoice);
  
  // this is for A1 and A2
  setup_Uint(cooUint1, apsi1, apsi2, NFac, dim, mu[0], intParam, intChoice,
	     SYS_TWO_COMPONENT, 1);
  setup_Uint(cooUint2, apsi1, apsi2, NFac, dim, mu[1], intParam, intChoice,
	     SYS_TWO_COMPONENT, 2);

  // this is for B1 and B2
  setup_Bint(cooB1, apsi1, apsi2, NFac, dim, intParam, intChoice,
	     SYS_TWO_COMPONENT, 1);
  setup_Bint(cooB2, apsi1, apsi2, NFac, dim, intParam, intChoice,
	     SYS_TWO_COMPONENT, 2);

  // this is for C
  setup_C(cooC, apsi1, apsi2, NFac, dim, intParam, intChoice);
  //  printf("setup terms.\n");
  
  coo_t *K, *M;
  K = coo_malloc(2*diff->nnz, 2*size, 2*size, flag);
  M = coo_malloc(2*diff->nnz+2*size, 2*size, 2*size, flag);

  setup_KM2C(K, diff, cooTrap, cooUint1, cooUint2, cooB1, cooB2,
	     cooC, 'K', diff->nnz, flag);
  setup_KM2C(M, diff, cooTrap, cooUint1, cooUint2, cooB1, cooB2,
	     cooC, 'M', diff->nnz, flag);

  //  printf("setup matrix.\n");

  // turn K and M from coo into rsb
  rsb_init(rsbTune);
  struct sp_lrep_t *LREP;
  LREP = splrep_malloc(N2, nev, subspaceSize, dim, K,
		       M, type, shift, maxIter, tol, flag);

  sp_setup_precond(LREP);

  //  printf("setup precond.\n");
  /*
  // tune and optimize thread count and matrix structure
  if(rsbTune == 1) {
    rsb_int_t tnK = 0;
    rsb_int_t tnM = 0;
    rsb_tune_SPMM(LREP->K, LREP->X, subspaceSize, 2*size, tnK, flag);
    rsb_tune_SPMM(LREP->M, LREP->X, subspaceSize, 2*size, tnM, flag);
    if(LREP->spPrecond->MDiag != NULL) {
      rsb_int_t tnDiag = 0;
      rsb_tune_SPMM(LREP->spPrecond->MDiag, LREP->X, subspaceSize,
		    size, tnDiag, flag);
    }
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    errval = rsb_lib_set_opt(RSB_IO_WANT_EXECUTING_THREADS, &tnK);
    if(errval != RSB_ERR_NO_ERROR) {
      rsb_perror(NULL, errval);
      printf("Tuning and Optimizing terminating with error.\n");
    }
  }
  */
  coo_free( M );
  coo_free( K );
  coo_free( cooC );
  coo_free( cooB1 );
  coo_free( cooB2 );
  coo_free( cooUint1 );
  coo_free( cooUint2 );
  coo_free( cooTrap );
  coo_free( diff );

  grid_free( grid );
  
  safe_free( apsi1 );
  safe_free( apsi2 );
  safe_free( wf1 );
  safe_free( wf2 );

  //  printf("setup LREP.\n");

  return LREP;
}
