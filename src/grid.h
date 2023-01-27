#pragma once

#include "gen.h"

// *******************************************************************
// Implements methods to build an euclidean grid, mostly for later use
// to build the differentiation operator and trap operator
//
// Includes the grid_t struct and following functions:
// - grid_malloc: allocates memory for grid_t
// - grid_free:   frees memory for grid_t
// - grid_setup:  sets up the grid in dim dimensions once it's allocated
// - idx_unravel: converts a flat index i into a couple of coordinate
//                arrays index
// *******************************************************************

typedef struct grid_t {
  int dim;
  int *N;
  double *ln;
  double **xn;
} grid_t;


struct grid_t* grid_malloc(const int dim, const int N[dim], const double ln[dim]);
void grid_free(grid_t *grid);
void grid_setup(grid_t *grid);
void* idx_unravel(int *index, int i, grid_t *grid);
