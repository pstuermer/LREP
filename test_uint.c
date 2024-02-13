#include "src/splrep.h"
#include "src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  const int dim = 3;
  int N[dim];

  N[0] = 4;
  N[1] = 4;
  N[2] = 4;

  // need to create some fake density
  // some fake mu
  // and some fake intParams

  double *apsi = (double*)xmalloc(N[0]*N[1]*N[2]*sizeof(double));
  for(int i = 0; i < N[0]*N[1]*N[2]; i++)
    apsi[i] = 1.0;

  double mu = 5.0 + 1.0e-12;

  double intParam[5] = { 0.0 };
  intParam[0] = 1.0;
  intParam[1] = 0.0;

  rsb_lib_init(RSB_NULL_INIT_OPTIONS);
  
  printf("Do bint test.\n");
  struct rsb_mtx_t *bint = NULL;
  bint = rsb_setup_int(apsi, 0.0, N, dim, intParam);
  
  printf("Do uint test.\n");
  struct rsb_mtx_t *uint = NULL;
  uint = rsb_setup_int(apsi, mu, N, dim, intParam);
  
  double *vals = xmalloc(N[0]*N[1]*N[2]*sizeof(double));
  int *row = xmalloc(N[0]*N[1]*N[2]*sizeof(int));
  int *col = xmalloc(N[0]*N[1]*N[2]*sizeof(int));

  rsb_mtx_get_coo(bint, vals, row, col, RSB_FLAG_C_INDICES_INTERFACE);

  print_array_d(N[0],1,vals);
  print_array(N[0],1, row);
  print_array(N[0],1, col);

  rsb_mtx_get_coo(uint, vals, row, col, RSB_FLAG_C_INDICES_INTERFACE);

  print_array_d(N[0]*N[1],1,vals);
  print_array(N[0]*N[1],1, row);
  print_array(N[0]*N[1],1, col);
  
  rsb_mtx_free( uint );
  rsb_mtx_free( bint );

  safe_free( apsi );
  safe_free( vals );
  safe_free( col );
  safe_free( row );
  
  rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
  return 0;
}
