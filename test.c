#include "src/splrep.h"
#include "src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  const double lx = 20.0;
  const int N = 8;

  double * diffMatrix = (double*)xmalloc(N*(N+1)/2*sizeof(double));

  //  print_array_d(N*N,1,diffMatrix);

  fourier_diff2(lx, N, diffMatrix);

  print_array_d(N*(N+1)/2, 1, diffMatrix);

  safe_free(diffMatrix);

  return 0;
}
