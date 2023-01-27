
#include "grid.h"

struct grid_t *grid_malloc(const int dim, const int N[dim], const double ln[dim]) {
  grid_t *grid = xmalloc(sizeof(grid_t));

  grid -> dim = dim;
  grid -> N = xmalloc(dim*sizeof(int));
  grid -> ln = xmalloc(dim*sizeof(double));
  grid -> xn = xmalloc(dim*sizeof(double *));

  for(int i = 0; i < dim; i++) {
    grid -> N[i] = N[i];
    grid -> ln[i] = ln[i];
    grid -> xn[i] = xmalloc(N[i]*sizeof(double));
  }

  return grid;
}

void grid_free(grid_t *grid) {
  assert(grid != NULL);

  for(int i = 0; i < grid->dim; i++) {
    safe_free( grid -> xn[i] );
  }

  safe_free( grid -> N );
  safe_free( grid -> ln );
  safe_free( grid -> xn );
  safe_free( grid );

  grid = NULL;
}

void grid_setup(grid_t *grid) {
  double dx;
  for(int i = 0; i < grid->dim; i++) {
    dx = (grid->ln[i])/(grid->N[i]);
    //    grid->N[i] -= 2;
    for(int j = 0; j < grid->N[i]; j++) {
      grid -> xn[i][j] = -(grid->ln[i])/2.0+(j+1)*dx;
    }
  }
}

void *idxunravel(int *index, int i, grid_t *grid) {
  if(grid->dim == 1) {
    index[0] = i;
  } else if (grid->dim == 2) {
    index[0] = i % (grid->N[0] - 1);
    index[1] = (i - index[0])/(grid->N[0]-1);
  } else if (grid->dim == 3) {
    index[0] = i % (grid->N[0] - 1);
    index[2] = i % ((grid->N[0])*(grid->N[1]-1));
    index[1] = (i - index[0] - index[2]*(grid->N[0])*(grid->N[1]))/(grid->N[0]-1);
  }
  return NULL;
}
