
#include "src/grid.h"
#include <rsb.h>

int main(void) {

  const int dim = 3;
  double lx[dim];
  int N[dim];

  grid_t *grid_1D = NULL, *grid_2D = NULL, *grid_3D = NULL;

  lx[0] = 2.0;
  lx[1] = 2.0;
  lx[2] = 2.0;
  N[0] = 16;
  N[1] = 16;
  N[2] = 64;

  grid_1D = grid_malloc(1, N, lx);
  grid_setup(grid_1D);
  print_array_d(N[0], 1, grid_1D->xn[0]);

  grid_2D = grid_malloc(2, N, lx);
  grid_setup(grid_2D);
  print_array_d(N[0], 1, grid_2D->xn[0]);
  print_array_d(N[1], 1, grid_2D->xn[1]);

  grid_3D = grid_malloc(dim, N, lx);
  grid_setup(grid_3D);
  print_array_d(N[0], 1, grid_3D->xn[0]);
  print_array_d(N[1], 1, grid_3D->xn[1]);
  print_array_d(N[2], 1, grid_3D->xn[2]);

  grid_free( grid_1D );
  grid_free( grid_2D );
  grid_free( grid_3D );
  
  return 0;
}
