#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <complex.h>
#include <time.h>
#include <string.h>

// *******************************************************************
// Includes a basic set of necessary libraries, as well as macros and
// utility functions that replace the common malloc and calloc with ones
// that have error checking included. PrintArray functions, print out
// an array to save myself from writing for-loops all the time.
// *******************************************************************

#ifndef M_PI
#define M_PI acos(-1.0)
#endif

#ifndef RAND_MAX
#define RAND_MAX((int) ((unsigned) ~0 >> 1))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef DATA_TYPE_DOUBLE
#define DATA_TYPE_DOUBLE 'D'
#endif

#ifndef DATA_TYPE_COMPLEX
#define DATA_TYPE_COMPLEX 'Z'
#endif

#ifndef SYS_ONE_COMPONENT
#define SYS_ONE_COMPONENT 1
#endif

#ifndef SYS_TWO_COMPONENT
#define SYS_TWO_COMPONENT 2
#endif

void *xmalloc(const int size);
void *xcalloc(const int num, const int size);
void safe_free(void *ptr);
void init_rand(void);
double real_rand(double low, double high);
double complex complex_rand(double complex low, double complex high);
void print_array_d(const int col, const int row, double *matrix);
void print_array_z(const int col, const int row, double complex *matrix);
void print_array(const int col, const int row, int *matrix);
