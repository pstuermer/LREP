
#include "rsbwrapper.h"
#include <rsb.h>
#include "openblaswrapper.h"

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
    mtxAp = rsb_mtx_alloc_from_coo_const(matrix->zval, matrix->row, matrix->col,
					 matrix->nnz, ztypeCode, matrix->rows,
					 matrix->cols, brA, bcA,
					 RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
					 | RSB_FLAG_DUPLICATES_SUM
					 | RSB_FLAG_DISCARD_ZEROS
					 | RSB_FLAG_HERMITIAN,
					 &errVal);
  if((!mtxAp) || (errVal != RSB_ERR_NO_ERROR)) {
    printf("Error allocating symmetric RSB matrix from coo: %d", errVal);
    rsb_perror(NULL, errVal);
  }

  return mtxAp;
}

struct rsb_mtx_t *rsb_mtx_from_coo_sym_new(double *vals, int *rows, int* cols,
					   int nnz, int ncols) {
  struct rsb_mtx_t *mtxAp = NULL;

  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t dtypeCode = RSB_NUMERICAL_TYPE_DEFAULT;

  mtxAp = rsb_mtx_alloc_from_coo_const(vals, rows, cols, nnz,
				       dtypeCode, ncols, ncols,
				       brA, bcA,
				       RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
				       | RSB_FLAG_DUPLICATES_SUM
				       | RSB_FLAG_DISCARD_ZEROS
				       | RSB_FLAG_SYMMETRIC,
				       &errVal);

  if((!mtxAp) || (errVal != RSB_ERR_NO_ERROR)) {
    printf("Error allocating symmetric RSB_matrix from arrays: %d\n", errVal);
    rsb_perror(NULL, errVal);
  }

  return mtxAp;
}


struct rsb_mtx_t *rsb_dadd(struct rsb_mtx_t *left,
			   struct rsb_mtx_t *right) {
  
  rsb_coo_idx_t left_nrow;
  rsb_coo_idx_t left_ncol;
  rsb_nnz_idx_t left_nnz;
  rsb_flags_t left_flag;
  rsb_type_t left_type;
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &left_nrow);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &left_ncol);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &left_nnz);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &left_flag);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &left_type);

  rsb_coo_idx_t right_nrow;
  rsb_coo_idx_t right_ncol;
  rsb_nnz_idx_t right_nnz;
  rsb_flags_t right_flag;
  rsb_type_t right_type;
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &right_nrow);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &right_ncol);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &right_nnz);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &right_flag);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &right_type);

  if(right_nrow != left_nrow)
    printf("left and right row count not equal in rsb_add_inplace.\n");
  if(right_ncol != left_ncol)
    printf("left and right col count not equal in rsb_add_inplace.\n");
  if(right_type != left_type)
    printf("trying to add matrices of different types in rsb_add_inplace.\n");

  int new_nnz = left_nnz + right_nnz;

  double *new_nz = xmalloc(new_nnz*sizeof(double));
  int *new_row = xmalloc(new_nnz*sizeof(int));
  int *new_col = xmalloc(new_nnz*sizeof(int));

  double *right_nz = xmalloc(right_nnz*sizeof(double));
  int *right_row = xmalloc(right_nnz*sizeof(int));
  int *right_col = xmalloc(right_nnz*sizeof(int));

  double *left_nz = xmalloc(left_nnz*sizeof(double));
  int *left_row = xmalloc(left_nnz*sizeof(int));
  int *left_col = xmalloc(left_nnz*sizeof(int));
  
  rsb_mtx_get_coo(left, left_nz, left_row, left_col, RSB_FLAG_C_INDICES_INTERFACE);
  rsb_mtx_get_coo(right, right_nz, right_row, right_col, RSB_FLAG_C_INDICES_INTERFACE);
  
#pragma omp parallel for schedule (static)
  for(int i = 0; i < left_nnz; i++) {
    new_nz[i] = left_nz[i];
    new_row[i] = left_row[i];
    new_col[i] = left_col[i];
  }

#pragma omp parallel for schedule (static)
  for(int i = 0; i < right_nnz; i++) {
    new_nz[i+left_nnz] = right_nz[i];
    new_row[i+left_nnz] = right_row[i];
    new_col[i+left_nnz] = right_col[i];
  }

  rsb_coo_sort(new_nz, new_row, new_col,
 	       left_nnz+right_nnz, right_nrow, right_ncol,
  	       RSB_NUMERICAL_TYPE_DOUBLE, RSB_FLAG_NOFLAGS);
  
  // allocate new matrix
  struct rsb_mtx_t *new = NULL;
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t dtypeCode = RSB_NUMERICAL_TYPE_DEFAULT;

  new = rsb_mtx_alloc_from_coo_const(new_nz, new_row, new_col,
				     new_nnz, dtypeCode, right_nrow,
				     right_ncol, brA, bcA,
				     RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
				     | RSB_FLAG_DUPLICATES_SUM
				     | RSB_FLAG_DISCARD_ZEROS
				     | RSB_FLAG_SYMMETRIC,
				     &errVal);

  if((!new) || (errVal != RSB_ERR_NO_ERROR)) {
      printf("Error allocating symmetric RSB matrix while inplace adding: %d", errVal);
      rsb_perror(NULL, errVal);
  }

  safe_free( new_nz );
  safe_free( new_row );
  safe_free( new_col );

  safe_free( right_nz );
  safe_free( right_row );
  safe_free( right_col );

  safe_free( left_nz );
  safe_free( left_row );
  safe_free( left_col );
  
  
  return new;
}

struct rsb_mtx_t *rsb_dsub(struct rsb_mtx_t *left,
			   struct rsb_mtx_t *right) {
  
  rsb_coo_idx_t left_nrow;
  rsb_coo_idx_t left_ncol;
  rsb_nnz_idx_t left_nnz;
  rsb_flags_t left_flag;
  rsb_type_t left_type;
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &left_nrow);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &left_ncol);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &left_nnz);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &left_flag);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &left_type);

  rsb_coo_idx_t right_nrow;
  rsb_coo_idx_t right_ncol;
  rsb_nnz_idx_t right_nnz;
  rsb_flags_t right_flag;
  rsb_type_t right_type;
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &right_nrow);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &right_ncol);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &right_nnz);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &right_flag);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &right_type);

  if(right_nrow != left_nrow)
    printf("left and right row count not equal in rsb_sub.\n");
  if(right_ncol != left_ncol)
    printf("left and right col count not equal in rsb_sub.\n");
  if(right_type != left_type)
    printf("trying to add matrices of different types in rsb_sub.\n");

  int new_nnz = left_nnz + right_nnz;

  double *new_nz = xmalloc(new_nnz*sizeof(double));
  int *new_row = xmalloc(new_nnz*sizeof(int));
  int *new_col = xmalloc(new_nnz*sizeof(int));

  double *right_nz = xmalloc(right_nnz*sizeof(double));
  int *right_row = xmalloc(right_nnz*sizeof(int));
  int *right_col = xmalloc(right_nnz*sizeof(int));

  double *left_nz = xmalloc(left_nnz*sizeof(double));
  int *left_row = xmalloc(left_nnz*sizeof(int));
  int *left_col = xmalloc(left_nnz*sizeof(int));

  rsb_mtx_get_coo(left, left_nz, left_row, left_col, RSB_FLAG_C_INDICES_INTERFACE);
  rsb_mtx_get_coo(right, right_nz, right_row, right_col, RSB_FLAG_C_INDICES_INTERFACE);

  // Parallelize this later
#pragma omp parallel for schedule(static)
  for(int i = 0; i < left_nnz; i++) {
    new_nz[i] = left_nz[i];
    new_row[i] = left_row[i];
    new_col[i] = left_col[i];
  }

#pragma omp parallel for schedule (static)
  for(int i = 0; i < right_nnz; i++) {
    new_nz[i+left_nnz] = right_nz[i]*(-1.0);
    new_row[i+left_nnz] = right_row[i];
    new_col[i+left_nnz] = right_col[i];
  }

  rsb_coo_sort(new_nz, new_row, new_col,
	       new_nnz, right_nrow, right_ncol,
	       RSB_NUMERICAL_TYPE_DOUBLE, RSB_FLAG_NOFLAGS);

  // allocate new matrix
  struct rsb_mtx_t *new = NULL;
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t dtypeCode = RSB_NUMERICAL_TYPE_DEFAULT;

  new = rsb_mtx_alloc_from_coo_const(new_nz, new_row, new_col,
				      new_nnz, dtypeCode, right_nrow,
				      right_ncol, brA, bcA,
				      RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
				      | RSB_FLAG_DUPLICATES_SUM
				      | RSB_FLAG_DISCARD_ZEROS
				      | RSB_FLAG_SYMMETRIC,
				      &errVal);

  if((!new) || (errVal != RSB_ERR_NO_ERROR)) {
      printf("Error allocating symmetric RSB matrix while inplace adding: %d", errVal);
      rsb_perror(NULL, errVal);
  }

  safe_free( new_nz );
  safe_free( new_row );
  safe_free( new_col );

  safe_free( right_nz );
  safe_free( right_row );
  safe_free( right_col );

  safe_free( left_nz );
  safe_free( left_row );
  safe_free( left_col );

  return new;
}

struct rsb_mtx_t *rsb_dkron(struct rsb_mtx_t *left,
			    struct rsb_mtx_t *right) {
  
  rsb_coo_idx_t left_nrow;
  rsb_coo_idx_t left_ncol;
  rsb_nnz_idx_t left_nnz;
  rsb_flags_t left_flag;
  rsb_type_t left_type;
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &left_nrow);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &left_ncol);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &left_nnz);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &left_flag);
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &left_type);

  rsb_coo_idx_t right_nrow;
  rsb_coo_idx_t right_ncol;
  rsb_nnz_idx_t right_nnz;
  rsb_flags_t right_flag;
  rsb_type_t right_type;
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T, &right_nrow);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T, &right_ncol);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &right_nnz);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T, &right_flag);
  rsb_mtx_get_info(right, RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T, &right_type);

  if(right_type != left_type)
    printf("trying to add matrices of different types in rsb_kron.\n");

  int new_nnz = left_nnz*right_nnz;
  int new_ncol = left_ncol*right_ncol;

  double *right_nz = xmalloc(right_nnz*sizeof(double));
  int *right_row = xmalloc(right_nnz*sizeof(int));
  int *right_col = xmalloc(right_nnz*sizeof(int));

  double *left_nz = xmalloc(left_nnz*sizeof(double));
  int *left_row = xmalloc(left_nnz*sizeof(int));
  int *left_col = xmalloc(left_nnz*sizeof(int));

  double *new_nz = xmalloc(new_nnz*sizeof(double));
  int *new_row = xmalloc(new_nnz*sizeof(int));
  int *new_col = xmalloc(new_nnz*sizeof(int));


  rsb_mtx_get_coo(left, left_nz, left_row, left_col, RSB_FLAG_C_INDICES_INTERFACE);
  rsb_mtx_get_coo(right, right_nz, right_row, right_col, RSB_FLAG_C_INDICES_INTERFACE);
  
#pragma omp parallel for schedule(static) collapse(2)
  for(int i = 0; i < left_nnz; i++) {
    for(int j = 0; j < right_nnz; j++) {
      new_nz[j+i*right_nnz] = left_nz[i] * right_nz[j];
      new_row[j+i*right_nnz] = right_nrow*left_row[i] + right_row[j];
      new_col[j+i*right_nnz] = right_ncol*left_col[i] + right_col[j];
    }
  }
  
  struct rsb_mtx_t *new = NULL;
  new = rsb_mtx_from_coo_sym_new(new_nz, new_row, new_col, new_nnz, new_ncol);

  /*  struct rsb_mtx_t *new = NULL;
  rsb_err_t errVal = RSB_ERR_NO_ERROR;
  new = rsb_mtx_alloc_from_coo_begin(new_nnz, RSB_NUMERICAL_TYPE_DOUBLE,
				     new_nrow, new_ncol,
				     RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS
				     | RSB_FLAG_DUPLICATES_SUM
				     | RSB_FLAG_DISCARD_ZEROS
				     | RSB_FLAG_SYMMETRIC,
				     &errVal);

  if((!new) || (errVal != RSB_ERR_NO_ERROR)) {
    printf("Error allocating new RSB matrix in dkron: %d", errVal);
    rsb_perror(NULL,errVal);
  }

  printf("started kronecker product.\n");
  for(int i = 0; i < left_nnz; i++) {

    double block_nz[right_nnz];
    int block_row[right_nnz];
    int block_col[right_nnz];

    int r = left_col[i];
    int s = left_row[i];

    // scale right_nz with current value of left
    cblas_dcopy(right_nnz, right_nz, 1, block_nz, 1);
    cblas_dscal(right_nnz, left_nz[i], block_nz, 1);

    // Parallelize this later
#pragma omp parallel for schedule(static)
    for(int j = 0; j < right_nnz; j++) {
      block_row[j] = right_nrow*s+right_row[j];
      block_col[j] = right_ncol*r+right_col[j];
    }

    rsb_mtx_set_vals(new, block_nz, block_row, block_col, right_nnz,
		     RSB_FLAG_C_INDICES_INTERFACE
		     | RSB_FLAG_DUPLICATES_KEEP_LAST);
  }

  rsb_mtx_alloc_from_coo_end(&new);
  */

  safe_free( new_nz );
  safe_free( new_row );
  safe_free( new_col );

  safe_free( right_nz );
  safe_free( right_row );
  safe_free( right_col );

  safe_free( left_nz );
  safe_free( left_row );
  safe_free( left_col );

  return new;
}

void rsb_SPMM(struct rsb_mtx_t *spMatrix, const void *dMatrix,
	      void *res, int nev, const char flag) {

  const RSB_DEFAULT_TYPE zero = 0;
  const RSB_DEFAULT_TYPE one = 1;
  const double complex zzero = 0;
  const double complex zone = 1;

  rsb_err_t errVal = RSB_ERR_NO_ERROR;

  switch(flag) {
  case 'Z':
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zzero,res,nev);
    break;

  default:
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zero,res,nev);
    break;
  }
  /*  
  if(flag == 'D') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zero,res,nev);

  } else if(flag == 'Z') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zzero,res,nev);
  }
  */
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

  switch(flag) {
  case 'Z':
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zminus,res,nev);
    break;

  default:
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &minus,res,nev);
    break;
  } 
  /*
  if(flag == 'D') {
     errVal = rsb_spmm(RSB_TRANSPOSITION_N,&one,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &minus,res,nev);
  } else if(flag == 'Z') {
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zone,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,
		      &zminus,res,nev);
		      }*/
  
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

  switch(flag) {
  case 'Z':
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &zbeta,res,nev);
    break;

  default:
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&dalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &dbeta,res,nev);
    break;
  }
  /*
  if(flag == 'D') {    
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&dalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &dbeta,res,nev);
  } else if(flag == 'Z') {  
    errVal = rsb_spmm(RSB_TRANSPOSITION_N,&zalpha,spMatrix,nev,
		      RSB_FLAG_WANT_ROW_MAJOR_ORDER,dMatrix,nev,    
		      &zbeta,res,nev);
		      }*/
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
