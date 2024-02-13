#include "src/splrep.h"
#include "src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  const int dim = 3;
  double lx[dim];
  int N[dim];

  lx[0] = 2.0;
  lx[1] = 2.0;
  lx[2] = 2.0;
  N[0] = 512;
  N[1] = 256;
  N[2] = 64;

  rsb_lib_init(RSB_NULL_INIT_OPTIONS);

  rsb_time_t dt;
  dt = -rsb_time();
  
  printf("Do 1D diff test.\n");
  struct rsb_mtx_t *diff1D = NULL;
  diff1D = rsb_setup_diff2(lx, N, 1);
  dt += rsb_time();
  printf("%lfs\n",dt);

  N[0] = 256;

  dt = -rsb_time();
  printf("Do 2D diff test.\n");
  struct rsb_mtx_t *diff2D = NULL;
  diff2D = rsb_setup_diff2(lx, N, 2);
  dt += rsb_time();
  printf("%lfs\n",dt);

  N[0] = 128;
  N[1] = 128;
  dt = -rsb_time();
  printf("Do 3D diff test.\n");
  struct rsb_mtx_t *diff3D = NULL;
  diff3D = rsb_setup_diff2(lx, N, 3);
  dt += rsb_time();
  printf("%lfs\n",dt);
  
  rsb_mtx_free( diff1D );
  rsb_mtx_free( diff2D );
  rsb_mtx_free( diff3D );

  rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
  return 0;
}
