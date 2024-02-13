#include "src/splrep.h"
#include "src/splopb4dcg.h"
#include "src/grid.h"
#include <rsb.h>

double omega = 1.0;
double V_trap_1D(double *pos) {
  return 0.5*omega*omega*pos[0]*pos[0];
}

double V_trap_2D(double *pos) {
  return 0.5*omega*omega*(pos[0]*pos[0]+
			  pos[1]*pos[1]);
}


double V_trap_3D(double *pos) {
  return 0.5*omega*omega*(pos[0]*pos[0]+
			  pos[1]*pos[1]+
			  pos[2]*pos[2]);
}

int main(void) {

  const int dim = 3;
  double lx[dim];
  int N[dim];

  lx[0] = 2.0;
  lx[1] = 2.0;
  lx[2] = 2.0;
  N[0] = 4;
  N[1] = 4;
  N[2] = 2;

  grid_t *grid_1D = NULL;
  grid_1D = grid_malloc(1, N, lx);
  grid_setup(grid_1D);

  grid_t *grid_2D = NULL;
  grid_2D = grid_malloc(2, N, lx);
  grid_setup(grid_2D);

  grid_t *grid_3D = NULL;
  grid_3D = grid_malloc(dim, N, lx);
  grid_setup(grid_3D);

  rsb_lib_init(RSB_NULL_INIT_OPTIONS);
  
  printf("Do trap 1D test.\n");
  struct rsb_mtx_t *trap_1D = NULL;
  trap_1D = rsb_setup_trap(V_trap_1D, grid_1D);

  printf("Do trap 2D test.\n");
  struct rsb_mtx_t *trap_2D = NULL;
  trap_2D = rsb_setup_trap(V_trap_2D, grid_2D);
  
  printf("Do trap 3D test.\n");
  struct rsb_mtx_t *trap_3D = NULL;
  trap_3D = rsb_setup_trap(V_trap_3D, grid_3D);
  
  double *vals = xmalloc(N[0]*N[1]*N[2]*sizeof(double));
  int *row = xmalloc(N[0]*N[1]*N[2]*sizeof(int));
  int *col = xmalloc(N[0]*N[1]*N[2]*sizeof(int));

  rsb_mtx_get_coo(trap_1D, vals, row, col, RSB_FLAG_C_INDICES_INTERFACE);

  print_array_d(N[0],1,vals);
  print_array(N[0],1, row);
  print_array(N[0],1, col);

  rsb_mtx_get_coo(trap_2D, vals, row, col, RSB_FLAG_C_INDICES_INTERFACE);
  
  print_array_d(N[0]*N[1],1,vals);
  print_array(N[0]*N[1],1, row);
  print_array(N[0]*N[1],1, col);
  
  rsb_mtx_get_coo(trap_3D, vals, row, col, RSB_FLAG_C_INDICES_INTERFACE);
  
  print_array_d(N[0]*N[1]*N[2],1,vals);
  print_array(N[0]*N[1]*N[2],1, row);
  print_array(N[0]*N[1]*N[2],1, col);
  
  rsb_mtx_free( trap_1D );
  rsb_mtx_free( trap_2D );
  rsb_mtx_free( trap_3D );
  grid_free( grid_1D );
  grid_free( grid_2D );
  grid_free( grid_3D );
  safe_free( vals );
  safe_free( row );
  safe_free( col );

  rsb_lib_exit(RSB_NULL_EXIT_OPTIONS);
  return 0;
}
