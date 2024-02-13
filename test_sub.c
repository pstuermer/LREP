#include "src/rsbwrapper.h"
#include "src/splrep.h"
#include "src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  struct rsb_mtx_t *left = NULL;
  const int bs = RSB_DEFAULT_BLOCKING;
  const int brA = bs, bcA = bs;
  rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
  rsb_err_t errval = RSB_ERR_NO_ERROR;
  const rsb_nnz_idx_t nnzA = 4;
  const rsb_coo_idx_t nrA = 3;
  const rsb_coo_idx_t ncA = 3;

  rsb_coo_idx_t IA[] = {0,1,2,2}; // nonzero row indices coordinates
  rsb_coo_idx_t JA[] = {0,1,2,2}; // nonzero colum indices coordinates
  RSB_DEFAULT_TYPE VA[] = {11,22,32,1}; // nonzero values


  struct rsb_mtx_t *right = NULL;
  const rsb_nnz_idx_t nnzB = 4;
  const rsb_coo_idx_t nrB = 3;
  const rsb_coo_idx_t ncB = 3;

  rsb_coo_idx_t IB[] = {0,1,2,2};
  rsb_coo_idx_t JB[] = {0,1,1,2};
  RSB_DEFAULT_TYPE VB[] = {1,2,3,4};

  printf("Hello, RSB!\n");
  printf("Initializing the library ...\n");

  if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) !=
     RSB_ERR_NO_ERROR)
    {
      printf("Error initializing the library!\n");
      //      goto err;
    }

  printf("Correctly initialized the library.\n");

  errval = RSB_ERR_NO_ERROR;

  left = rsb_mtx_alloc_from_coo_const(VA,IA,JA,nnzA,typecode,
				       nrA,ncA,brA,bcA,
				       RSB_FLAG_NOFLAGS
				       | RSB_FLAG_DUPLICATES_SUM,
				       &errval);
  if((!left) || (errval != RSB_ERR_NO_ERROR))
    {
      printf("Error while allocating the matrix!\n");
      //      goto err;
    }
  printf("Correctly allocated a matrix.\n");

  errval = RSB_ERR_NO_ERROR;

  right = rsb_mtx_alloc_from_coo_const(VB,IB,JB,nnzB,typecode,
				       nrB,ncB,brA,bcA,
				       RSB_FLAG_NOFLAGS
				       | RSB_FLAG_DUPLICATES_SUM,
				       &errval);
  if((!right) || (errval != RSB_ERR_NO_ERROR))
    {
      printf("Error while allocating the matrix!\n");
      //      goto err;
    }
  printf("Correctly allocated a matrix.\n");

  left = rsb_dsub(left, right);

  rsb_nnz_idx_t new_nnz;
  rsb_mtx_get_info(left, RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T, &new_nnz);

  rsb_coo_idx_t IN[new_nnz];
  rsb_coo_idx_t JN[new_nnz];
  RSB_DEFAULT_TYPE VN[new_nnz];

  rsb_mtx_get_coo(left, VN, IN, JN, RSB_FLAG_C_INDICES_INTERFACE);

  print_array(new_nnz, 1, IN);
  print_array(new_nnz, 1, JN);
  print_array_d(new_nnz, 1, VN);
  rsb_mtx_free(left);
  rsb_mtx_free(right);


  return 0;

  /* err:
  rsb_perror(NULL,errval);
  printf("Program terminating with error.\n");
  return -1;*/
}
