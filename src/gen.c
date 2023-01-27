#include "gen.h"

void *xmalloc(const int size) {
  void* ptr = malloc(size);
  if (ptr == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }

  return ptr;
}

void *xcalloc(const int num, const int size) {
  void *ptr = calloc(num, size);
  if(ptr == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }

  return ptr;
}

void safe_free(void *ptr) {
  if(ptr != NULL) {

    free( ptr );
    ptr = NULL;
  }
}

void init_rand(void) {
  srand((int) time(NULL));
}

double real_rand(double low, double high) {
  double d;
  d = (double) rand() / ((double) RAND_MAX + 1);
  return (low + d * (high - low));
}

double complex complex_rand(double complex low, double complex high) {
  double complex d;
  d = (double complex) rand() / ((double complex) RAND_MAX + 1);
  return (low + d * (high - low));
}

void print_array_d(const int col, const int row, double *matrix) {
  for(int i = 0; i < row; i++) {
    for(int j = 0; j < col; j++) {
      printf("%.5e\t", matrix[j+i*col]);
    }
    printf("\n");
  }
  printf("\n");
}

void print_array_z(const int col, const int row, double complex *matrix) {
  for(int i = 0; i < row; i++) {
    for(int j = 0; j < col; j++) {
      printf("%.5e + %.5e i\t", creal(matrix[j+i*col]), cimag(matrix[j+i*col]));
    }
    printf("\n");
  }
  printf("\n");
}

void print_array(const int col, const int row, int *matrix) {
  for(int i = 0; i < row; i++) {
    for(int j = 0; j < col; j++) {
      printf("%d\t", matrix[j+i*col]);
    }
    printf("\n");
  }
  printf("\n");
}

