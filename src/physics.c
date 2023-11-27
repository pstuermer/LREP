
#include "physics.h"
#include <omp.h>

struct coo_t *setup_diff2(const double *ln, const int *N, const int dim,
			 const char flag) {
  coo_t *diff;
  if(dim == 1) {
    diff = diff_1D(ln, N, flag);
    return diff;
  } else if (dim == 2) {
    diff = diff_2D(ln, N, flag);
    return diff;
  } else if (dim == 3) {
    diff = diff_3D(ln, N, flag);
    return diff;
  } else {
    perror(0);
    return 0;
  }
}

struct coo_t *diff_1D(const double *ln, const int *N, const char flag) {
  double *diff;
  diff = xmalloc((N[0]*(N[0]+1)/2.0) * sizeof(double));
  //  const int N1 = N[0]+1;
  fourier_diff2(ln[0], N[0], diff);
  
  coo_t* diffFin;
  diffFin = coo_malloc((N[0]*(N[0]+1)/2.0), N[0], N[0], flag);

  int index = 0;
  
  for(int i = 0; i < N[0]; i++) {
    for(int j = i; j < N[0]; j++) {
      index = j+i*N[0] - (i*(i+1))/2;
      if(diff[index] == 0)
	continue;
      else
	coo_insert(diffFin, diff[index], 0, j, i);
    }
  }
  
  safe_free( diff );
  
  return diffFin;
}

struct coo_t *diff_2D(const double *ln, const int *N, const char flag) { 
  
  double *diff1, *diff2;
  diff1 = xmalloc(N[0]*N[0]*sizeof(double));
  diff2 = xmalloc(N[1]*N[1]*sizeof(double));

  //  int N1[2];
  //  for(int i = 0;i < 2;i++)
  //    N1[i] = N[i]+1;

  fourier_diff2(ln[0], N[0], diff1);
  fourier_diff2(ln[1], N[1], diff2);
  
  coo_t *kron1, *kron2, *unit1, *unit2, *diffFin1, *diffFin2;
  diffFin1 = coo_malloc(N[0]*N[0], N[0], N[0], flag);
  diffFin2 = coo_malloc(N[1]*N[1], N[1], N[1], flag);
  unit1 = coo_malloc(N[0], N[0], N[0], flag);
  unit2 = coo_malloc(N[1], N[1], N[1], flag);
  kron1 = coo_malloc(diffFin1->nnz*unit1->nnz,   
		     diffFin1->rows*unit1->rows,
		     diffFin1->cols*unit1->cols, flag);
  kron2 = coo_malloc(diffFin2->nnz*unit2->nnz,
		     diffFin2->rows*unit2->rows,
		     diffFin2->cols*unit2->cols, flag);

  // also need to allocate everything here 
  for(int i = 0; i < N[0]; i++)
    coo_insert(unit1, 1.0, 0, i, i);
  for(int i = 0; i < N[1]; i++) 
    coo_insert(unit2, 1.0, 0, i, i);

  for(int i = 0;i < N[0];i++) {  
    for(int j = 0;j < N[0];j++) { 
      if (diff1[j+i*N[0]] == 0)   
	continue;
      else   
	coo_insert(diffFin1, diff1[j+i*N[0]], 0, j, i);
    }  
  }
  for(int i = 0;i < N[1];i++) {  
    for(int j = 0;j < N[1];j++) {
      if (diff2[j+i*N[0]] == 0)   
	continue;
      else 
	coo_insert(diffFin2, diff2[j+i*N[1]], 0, j, i);
    }
  }  
  coo_kron(unit1, diffFin2, kron1);
  coo_kron(diffFin1, unit2, kron2);

  coo_t *diffFin;
  diffFin = coo_malloc(kron1->nnz+kron2->nnz,   
		       kron1->rows, kron1->cols,flag);
  
  coo_add(kron1, kron2, diffFin);
  
  safe_free( diff1 );
  safe_free( diff2 );
  coo_free( diffFin1 );
  coo_free( diffFin2 );
  coo_free( unit1 );
  coo_free( unit2 );
  coo_free( kron1 );
  coo_free( kron2 );
  return diffFin;
}


struct coo_t *diff_3D(const double *ln, const int *N, const char flag) {
  double *diff1, *diff2, *diff3;
  diff1 = xmalloc(N[0]*N[0]*sizeof(double));
  diff2 = xmalloc(N[1]*N[1]*sizeof(double));
  diff3 = xmalloc(N[2]*N[2]*sizeof(double));
  //  int N1[3];
  //  for(int i = 0;i < 3;i++)   
  //    N1[i] = N[i]+1;
  
  fourier_diff2(ln[0], N[0], diff1);
  fourier_diff2(ln[1], N[1], diff2);
  fourier_diff2(ln[2], N[2], diff3);
  
  coo_t *kron1, *kron2, *kron3, *unit1, *unit2, *unit3,   
    *diffFin1, *diffFin2, *diffFin3, *work1, *work2, *work3;
  
  diffFin1 = coo_malloc(N[0]*N[0], N[0], N[0], flag);
  diffFin2 = coo_malloc(N[1]*N[1], N[1], N[1], flag);
  diffFin3 = coo_malloc(N[2]*N[2], N[2], N[2], flag);
  unit1 = coo_malloc(N[0], N[0], N[0], flag);
  unit2 = coo_malloc(N[1], N[1], N[1], flag);
  unit3 = coo_malloc(N[2], N[2], N[2], flag);
  
  work1 = coo_malloc(diffFin1->nnz*unit2->nnz,
		     diffFin1->rows*unit2->rows,
		     diffFin1->cols*unit2->cols, flag);
  work2 = coo_malloc(diffFin2->nnz*unit1->nnz,
		     diffFin2->rows*unit1->rows,
		     diffFin2->cols*unit1->cols, flag);
  work3 = coo_malloc(diffFin3->nnz*unit2->nnz,
		     diffFin3->rows*unit2->rows,   
		     diffFin3->cols*unit2->cols, flag);
  
  kron1 = coo_malloc(work1->nnz*unit3->nnz,
		     work1->rows*unit3->rows,
		     work1->cols*unit3->cols, flag);
  kron2 = coo_malloc(work2->nnz*unit3->nnz,
		     work2->rows*unit3->rows,
		     work2->cols*unit3->cols, flag);
  kron3 = coo_malloc(work3->nnz*unit1->nnz,
		     work3->rows*unit1->rows,   
		     work3->cols*unit1->cols, flag);
  
  // allocate everything
  for(int i = 0;i < N[0];i++) 
    coo_insert(unit1, 1.0, 0, i, i);
  for(int i = 0;i < N[1];i++)
    coo_insert(unit2, 1.0, 0, i, i);
  for(int i = 0;i < N[2];i++) 
    coo_insert(unit3, 1.0, 0, i, i);
  
  for(int i = 0;i < N[0];i++) {   
    for(int j = 0;j < N[0];j++) { 
      if (diff1[j+i*N[0]] == 0)
	continue;
      else  
	coo_insert(diffFin1, diff1[j+i*N[0]], 0, j, i);
    }   
  } 
  
  for(int i = 0;i < N[1];i++) {   
    for(int j = 0;j < N[1];j++) { 
      if (diff2[j+i*N[1]] == 0)
	continue;
      else  
	coo_insert(diffFin2, diff2[j+i*N[1]], 0, j, i);
    }   
  } 
  
  for(int i = 0;i < N[2];i++) {   
    for(int j = 0;j < N[2];j++) { 
      if (diff3[j+i*N[2]] == 0)
	continue;
      else  
	coo_insert(diffFin3, diff2[j+i*N[2]], 0, j, i);
    }   
  } 
  
  coo_kron(diffFin1, unit2, work1);
  coo_kron(work1, unit3, kron1);
  
  coo_kron(unit1, diffFin2, work2);
  coo_kron(work2, unit3, kron2);
  
  coo_kron(unit2, diffFin3, work3);
  coo_kron(unit2, work3, kron3);
  
  coo_t *diffWork, *diffFin;
  diffWork = coo_malloc(kron1->nnz+kron2->nnz,   
			kron1->rows, kron1->cols, flag);
  diffFin = coo_malloc(diffWork->nnz+kron3->nnz,   
		       kron3->rows, kron3->cols, flag);
  coo_add(kron1, kron2, diffWork);
  coo_add(diffWork, kron3, diffFin);
  
  coo_free( diffWork );
  coo_free( kron3 );
  coo_free( kron2 );
  coo_free( kron1 );
  coo_free( work3 );
  coo_free( work2 );
  coo_free( work1 );
  coo_free( diffFin3 );
  coo_free( diffFin2 );
  coo_free( diffFin1 );
  coo_free( unit3 );
  coo_free( unit2 );
  coo_free( unit1 );
  safe_free( diff1 );
  safe_free( diff2 );
  safe_free( diff3 );
  
  return diffFin;
}

void fourier_diff2(const double lx, const int N, double *diffMatrix) { 
  double h, work;
  int l = 0;
  int index = 0;
  h = lx/N;

#pragma omp parallel for private(work) shared(diffMatrix) schedule (static)  
  // loop through matrix elements   
  for (int n = 1; n < N+1; n++) {   
    for (int j = n; j < N+1; j++) {
      index = (j-1)+(n-1)*N - (n*(n-1))/2;
      work = 0.0;
      diffMatrix[index] = 0.0;
      // sum up different sin components. There is probably 
      // a faster way to do this
      for (int k = 1; k < N+1; k++) {
	l = k - N/2;
	work = l*l*sin(2*M_PI*l*j/N)*sin(2*M_PI*l*n/N);
	work += l*l*cos(2*M_PI*l*j/N)*cos(2*M_PI*l*n/N);
	
	diffMatrix[index] += (-8.0*M_PI*M_PI)/(N*N*N*h*h)*work;
      }
      diffMatrix[index] *= -0.25;
    }   
  } 
}

void fourier_diff1(const double lx, const int N, double *diffMatrix) {
  double h, work;
  int l = 0;
  h = lx/N;

  // loop through matrix elements
  for (int n = 1; n < N+1; n++) {
    for (int j = 1; j < N+1; j++) {
      work = 0.0;
      diffMatrix[(j-1)+(n-1)*N] = 0.0;
      // sum up different fourier components. There is probably
      // a faster way to do this
      for (int k = 0; k < N+1; k++) {
	l = k-N/2;
	work = l*sin(2*M_PI*l*j/N)*cos(2*M_PI*l*n/N);
	work -= l*cos(2*M_PI*l*j/N)*sin(2*M_PI*l*n/N);
	diffMatrix[(j-1)+(n-1)*N] += 2*M_PI/(N*N*h)*work;
      }
      diffMatrix[(j-1)+(n-1)*N] *= 0.25;
    }
  }
}
	

void load_wf(double complex *wf, const int N, 
	     const int factor, char fileName[]) {   
  double complex z = 0.0 + I*0.0;
  double buffer = 0.0;
  int counter = 0;
  int isReal = 1;

  double complex *tempWf;
  tempWf = xmalloc(N*sizeof(double complex));
  
  FILE* inStream = fopen(fileName, "r");
  if (inStream) {   
    while (fscanf(inStream, "%lf", &buffer) && counter < N) {   
      if (isReal) { 
	isReal = 0;
	z = CMPLX(buffer, cimag(z));
      } else {  
	isReal = 1;
	z = CMPLX(creal(z), buffer);
	tempWf[counter] = z;
	counter++;
      } 
    }   
    fclose(inStream);
  } else {  
    /* Provides some error diagnostic */
    fprintf(stderr, "Could not open %s: ", fileName);
    perror(0);
    errno = 0;
  } 

  counter = 0;
  for(int i = 0; i < N; i++) {
    if(!(N%factor)) {
      wf[counter] = tempWf[i];
      counter++;
    }
  }

  safe_free( tempWf );
}

void load_wf_2C(double complex *wf1, double complex *wf2,
		const int N, const int factor, char fileName[]) {
  double complex z = 0.0 + I*0.0;
  double buffer = 0.0;
  int counter = 0;
  int wfCounter = 1;
  int isReal = 1;

  double complex *tempWf1, *tempWf2;
  tempWf1 = xmalloc(N*sizeof(double complex));
  tempWf2 = xmalloc(N*sizeof(double complex));

  FILE *inStream = fopen(fileName, "r");
  if (inStream) {
    while (fscanf(inStream, "%lf", &buffer) && counter < N) {
      if (isReal && wfCounter) {
	isReal = 0;
	z = CMPLX(buffer, cimag(z));
      } else if (!isReal && wfCounter) {
	isReal = 1;
	wfCounter = 0;
	z = CMPLX(creal(z), buffer);
	tempWf1[counter] = z;
      } else if (isReal && !wfCounter) {
	isReal = 0;
	z = CMPLX(buffer, cimag(z));
      } else if (!isReal && !wfCounter) {
	isReal = 1;
	wfCounter = 1;
	z = CMPLX(creal(z), buffer);
	tempWf2[counter] = z;
	counter++;
      }
    }
    fclose(inStream);
  } else {
    // Provides some error diagnostic
    fprintf(stderr, "Could not open %s: ", fileName);
    perror(0);
    errno = 0;
  }

  counter = 0;
  for(int i = 0; i < N; i++) {
    if(!(N%factor)) {
      wf1[counter] = tempWf1[i];
      wf2[counter] = tempWf2[i];
      counter++;
    }
  }

  safe_free( tempWf1 );
  safe_free( tempWf2 );


}
  
double *get_apsi(double complex *wf, const int N) { 
  double *apsi;
  apsi = xmalloc(N*sizeof(double));
  for(int i = 0;i < N;i++) {  
    apsi[i] = sqrt(creal(wf[i])*creal(wf[i]) + cimag(wf[i])*cimag(wf[i]));
    if(apsi[i] < 1.0e-15)
      apsi[i] = 0.0;
  } 
  return apsi;
}

double complex *get_psi_sq(double complex *wf, const int N) {
  double complex *psiSq;
  psiSq  = xmalloc(N*sizeof(double complex));
  for(int i = 0; i < N; i++) {
    psiSq[i] = wf[i];
  }
  return psiSq;
}

double complex *get_ccpsi_sq(double complex *wf, const int N) {
  double complex *ccPsiSq;
  ccPsiSq = xmalloc(N*sizeof(double complex));
  for(int i = 0; i < N; i++) {
    ccPsiSq[i] = conj(wf[i]);
  }
  return ccPsiSq;
}

double *cut_apsi(double *apsi, const int *N, const int dim) {   
  double *newApsi;
  int counter = 0;
  if(dim == 1) {
    newApsi = xmalloc((N[0]-2)*sizeof(double));
    for(int i = 1;i < N[0]-1;i++) {   
      newApsi[counter] = apsi[i];
      counter++;
    }
    return newApsi;
  } else if (dim == 2) {
    newApsi = xmalloc((N[0]-2)*(N[1]-2)*sizeof(double));
    for(int i = 1;i < N[0]-1;i++) {   
      for(int j = 1;j < N[1]-1;j++) { 
	newApsi[counter] = apsi[j+ i*N[1]];
	counter++;
      } 
    }  
    return newApsi;
  } else if (dim == 3) {
    newApsi = xmalloc((N[0]-2)*(N[1]-2)*(N[2]-2)*sizeof(double));
    for(int i = 1;i < N[0]-1;i++) {   
      for(int j = 1;j < N[1]-1;j++) { 
	for(int k = 1;k < N[2]-1;k++) {   
	  newApsi[counter] = apsi[k + j*(N[1]) + i*(N[1])*(N[2])];
	  counter++;
	}   
      } 
    }  
    return newApsi;
  }
  
  return 0;
}

double complex *cut_wf(double complex *wf, const int *N, const int dim) {
  double complex *newWf;
  int counter = 0;
  if(dim == 1) {
    newWf = xmalloc((N[0]-2)*sizeof(double complex));
    for(int i = 1; i < N[0]-1; i++) {
      newWf[counter] = wf[i];
      counter++;
    }
    return newWf;
  } else if(dim == 2) {
    newWf = xmalloc((N[0]-2)*(N[1]-2)*sizeof(double complex));
    for(int i = 1; i < N[0]-1; i++) {
      for(int j = 1; j < N[1]-1; j++) {
	newWf[counter] = wf[j+i*N[1]];
	counter++;
      }
    }
    return newWf;
  } else if(dim == 3) {
    newWf = xmalloc((N[0]-2)*(N[1]-2)*(N[2]-2)*sizeof(double complex));
    for(int i = 1; i < N[0]-1; i++) {
      for(int j = 1; j < N[1]-1; j++) {
	for(int k = 1; k < N[2]-1; k++) {
	  newWf[counter] = wf[k+j*N[1]+i*N[1]*N[2]];
	  counter++;
	}
      }
    }
    return newWf;
  }
  return 0;
}

/*void setup_mu(double *diffMatrix, double* apsi1, double* apsi2, double mu[2],
	      const int intChoice, const double *intParam, grid_t sys) {

  // get wf
  int size = 1;
  for (i = 0; i < grid_t->dim; i++)
    size *= grid_t->N[i];
  
  double *wf1, *wf2;
  wf1 = xmalloc(

  // get gradient of wf

  // get mu
  */

void setup_trap(coo_t *cooTrap, const int *N, const double *param, grid_t *grid,
		const int dim, const int trapChoice) { 
  int index[dim];
  
  if(dim == 1) {
    for(int i = 0;i < N[0];i++) { 
      index[0] = i;
      coo_insert(cooTrap, trap_1D(param, grid, index, trapChoice),
		 0, i, i);
    }   
    
  } else if (dim == 2) {
    for(int i = 0;i < N[0];i++) { 
      for(int j = 0;j < N[1];j++) {   
	index[0] = i;
	index[1] = j;
	coo_insert(cooTrap, trap_2D(param, grid, index, trapChoice), 
		   0, j+i*N[1], j+i*N[1]);
      } 
    }   
    
  } else if (dim == 3) {
    for(int i = 0;i < N[0];i++) { 
      for(int j = 0;j < N[1];j++) {   
	for(int k = 0;k < N[2];k++) { 
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  coo_insert(cooTrap, trap_3D(param, grid, index, trapChoice),   
		     0, k+j*N[2]+i*N[2]*N[1], k+j*N[2]+i*N[2]*N[1]);
	}   
      } 
    }   
    
  } else {  
    printf("Dimensionality error in V_trap.\n");
    perror(0);
  } 

}   

double trap_1D(const double *param, grid_t *sys, int *index, const int trapChoice) {  
  
  double x;
  x = sys->xn[0][index[0]];
  
  if(!trapChoice) { 
    return 0.0;// do nothing as no trap chosen 
  } else if (trapChoice == 1) { 
    
    return 0.5*param[0]*param[0]*x*x;
  } else {  
    printf("Invalid trapchoice.\n");
    perror(0);
  } 
  
  return 0.0;
}   

double trap_2D(const double *param, grid_t *sys, int *index, const int trapChoice) {  
  
  double x, y;
  x = sys->xn[0][index[0]];
  y = sys->xn[1][index[1]];
  
  if(!trapChoice) { 
    return 0.0;// do nothing as no trap chosen 
  } else if (trapChoice == 1) { 
    
    return 0.5*param[0]*param[0]*(x*x+y*y);
  } else {  
    printf("Invalid trapchoice.\n");
    perror(0);
  } 
  
  return 0.0;
}

double trap_3D(const double *param, grid_t *sys, int *index, const int trapChoice) {  
  
  double x, y, z;
  x = sys->xn[0][index[0]];
  y = sys->xn[1][index[1]];
  z = sys->xn[2][index[2]];
  
  if(!trapChoice) { 
    return 0.0;// do nothing as no trap chosen 
  } else if (trapChoice == 1) { 
    
    return 0.5*param[0]*param[0]*(x*x+y*y+z*z);
  } else {  
    printf("Invalid trapchoice.\n");
    perror(0);
  } 
  
  return 0.0;
}   

void setup_Uint(coo_t *cooUint, double *apsi1, double *apsi2, const int *N,
		const int dim, const double mu, const double *intParam,
		const int intChoice, const int sysComp, int comp) {
  if(sysComp == SYS_ONE_COMPONENT) {
    setup_Uint1C(cooUint, apsi1, N, dim, mu, intParam, intChoice);
  } else if (sysComp == SYS_TWO_COMPONENT) {
    setup_Uint2C(cooUint, apsi1, apsi2, N, dim, mu, intParam, intChoice, comp);
  }

}

void setup_Uint1C(coo_t *cooUint, double *apsi, const int *N, const int dim,
		  const double mu, const double *intParam, const int intChoice){
  double dapsi;
  int size = 1;
   for(int i = 0; i < dim; i++)
    size *= N[i];
  for(int i = 0; i < size; i++) {
    dapsi = apsi[i];
    coo_insert(cooUint, d_Uint1C(dapsi, mu, intParam, intChoice), 0.0, i, i);
  }

}

void setup_Uint2C(coo_t *cooUint, double *apsi1, double *apsi2, const int *N,
		  const int dim, const double mu, const double *intParam,
		  const int intChoice, int comp) {
  int size = 1;
  double dapsi1, dapsi2;
  for(int i = 0; i < dim; i++) 
    size *= N[i];
  for(int i = 0; i < size; i++) {
    dapsi1 = apsi1[i];
    dapsi2 = apsi2[i];
    coo_insert(cooUint, d_Uint2C(dapsi1, dapsi2, mu, intParam, intChoice, comp),
	       0.0, i, i);
  }

}


double d_Uint1C(double apsi, const double mu, const double *intParam, const int intChoice) {   
  if(!intChoice) {  
    return -mu;
  } else if (intChoice == 1) {  
    // MF   
    return 2.0*intParam[0]*apsi*apsi-mu;
  } else if (intChoice == 2) {  
    // 1D-LHY
    return 2.0*intParam[0]*apsi*apsi+1.5*intParam[1]*apsi-mu;
  } else if (intChoice == 3) {
    // 2D-LHY
    return apsi*apsi*(2.0*log(apsi*apsi)+1.0)-mu;
  } else {  
    printf("Invalid Uintchoice.\n");
    perror(0);
  } 
  
  return 0.0;
}

double d_Uint2C(double apsi1, double apsi2, const double mu,
		const double *intParam, const int intChoice, int comp) {
  double n1 = apsi1*apsi1;
  double n2 = apsi2*apsi2;
  if(!intChoice) {
    return -mu;
  } else if (intChoice == 1 && comp == 1) {
    // MF
    return 2*intParam[0]*apsi1*apsi1 + intParam[2]*apsi2*apsi2 - mu;
  } else if (intChoice == 1 && comp == 2) {
    // MF
    return 2*intParam[1]*apsi2*apsi2 + intParam[2]*apsi2*apsi2 - mu;
    

  } else if (intChoice == 3 && comp == 1) {
    return 2.0*intParam[1]*n1/intParam[0] +
      intParam[2]*n2 - intParam[3]*intParam[1]*sqrt(intParam[1])/(M_PI*intParam[0])*
      (3.0*n1/intParam[0]+2.0*n2*intParam[0])/
      sqrt(4.0*(n1/intParam[0]+n2*intParam[0])) - mu;
  } else if (intChoice == 3 && comp == 2) {
    return 2.0*intParam[1]*n2*intParam[0] +
      intParam[2]*n1 - intParam[3]*intParam[1]*sqrt(intParam[1])*intParam[0]/M_PI*
      (2.0*n1/intParam[0]+3.0*n2*intParam[0])/
      sqrt(4.0*(n1/intParam[0]+n2*intParam[0])) - mu;

    
  } else if (intChoice == 4 && comp == 1) {
    return (2.0*intParam[0]-3.0*intParam[0]*sqrt(intParam[0])/
	    (2.0*M_PI*sqrt(n1+n2)))*n1+
      (intParam[1]-intParam[0]-2.0*intParam[0]*sqrt(intParam[0])/
       (2.0*M_PI*sqrt(n1+n2)))*n2-mu;
  } else if (intChoice == 4 && comp == 2) {
    return (2.0*intParam[0]-3.0*intParam[0]*sqrt(intParam[0])/
	    (2.0*M_PI*sqrt(n1+n2)))*n2+
      (intParam[1]-intParam[0]-2.0*intParam[0]*sqrt(intParam[0])/
       (2.0*M_PI*sqrt(n1+n2)))*n1-mu;
    
  } else {
    printf("Invalid UintChoice.\n");
    perror(0);
  }
  return 0.0;
}

void setup_Bint(coo_t *coo_Bint, double *apsi1, double *apsi2,
		const int *N, const int dim, const double *intParam,
		const int intChoice, const int sysComp, int comp) {
  if(sysComp == SYS_ONE_COMPONENT)
    setup_Bint1C(coo_Bint, apsi1, N, dim, intParam, intChoice);
  else if(sysComp == SYS_TWO_COMPONENT)
    setup_Bint2C(coo_Bint, apsi1, apsi2, N, dim, intParam, intChoice, comp);

}

void setup_Bint1C(coo_t *coo_Bint, double *apsi, const int *N,
		  const int dim, const double *intParam, const int intChoice) {
  int size = 1;
  double dapsi;
  
  for(int i = 0; i < dim; i++)
    size *= N[i];

  for(int i = 0; i < size; i++) {
    dapsi = apsi[i];
    coo_insert(coo_Bint, d_Bint1C(dapsi, intParam, intChoice), 0.0, i, i);
  }

}

void setup_Bint2C(coo_t *coo_Bint, double* apsi1, double *apsi2,
		  const int *N, const int dim, const double *intParam, 
		  const int intChoice, int comp) {
  int size = 1;
  double dapsi1, dapsi2;
  for(int i = 0; i < dim; i++)
    size *= N[i];

  for(int i = 0; i < size; i++) {
    dapsi1 = apsi1[i];
    dapsi2 = apsi2[i];
    coo_insert(coo_Bint, d_Bint2C(dapsi1, dapsi2, intParam, intChoice, comp),
	       0.0, i, i);
  }

}

/*
void *zsetup_Bint(coo_t *coo_Bint1, double complex *psiSq,
		  double complex *ccPsiSq, const int *N, const int dim, 
		  const double *intParam, const int intChoice) {

  int index;
  double complex zpsiSq;//, zccPsiSq;

  if(dim == 1) {
    for(int i = 0; i < N[0]; i++) {
      zpsiSq = ccPsiSq[i];
      zpsiSq = psiSq[i];
      //      zccPsiSq = ccPsiSq[i];
      coo_insert(coo_Bint1, creal(z_Bint(zpsiSq, intParam, intChoice)),
		 cimag(z_Bint(zpsiSq, intParam, intChoice)), i, i);
      //      coo_insert(coo_Bint2, creal(z_Bint(zccPsiSq, intParam, intChoice)),
      //		 cimag(z_Bint(zccPsiSq, intParam, intChoice)), i, i);
    }
  } else if(dim == 2) {
    for(int i = 0; i < N[0]; i++) {
      for(int j = 0; j < N[1]; j++) {
	index = j+i*N[1];
	zpsiSq = ccPsiSq[index];
	zpsiSq = psiSq[index];
	//	zccPsiSq = ccPsiSq[i];
	coo_insert(coo_Bint1, creal(z_Bint(zpsiSq, intParam, intChoice)),
		   cimag(z_Bint(zpsiSq, intParam, intChoice)), index, index);
	//	coo_insert(coo_Bint2, creal(z_Bint(zccPsiSq, intParam, intChoice)),
	//		   cimag(z_Bint(zccPsiSq, intParam, intChoice)), index, index);
      }
    }
  } else if(dim == 3) {
    for(int i = 0; i < N[0]; i++) {
      for(int j = 0; j < N[1]; j++) {
	for(int k = 0; k < N[2]; k++) {
	  index = k+j*N[2]+i*N[1]*N[2];
	  zpsiSq = ccPsiSq[index];
	  zpsiSq = psiSq[index];
	  //	  zccPsiSq = ccPsiSq[index];
	  coo_insert(coo_Bint1, creal(z_Bint(zpsiSq, intParam, intChoice)),
		     cimag(z_Bint(zpsiSq, intParam, intChoice)), index, index);
	  //	  coo_insert(coo_Bint2, creal(z_Bint(zccPsiSq, intParam, intChoice)),
	  //		     cimag(z_Bint(zccPsiSq, intParam, intChoice)), index, index);
	}
      }
    }
  }

  return NULL;
}
*/
double d_Bint1C(double apsi, const double *intParam, const int intChoice) { 
    if(!intChoice) {  
      return 0.0;
    } else if (intChoice == 1) {  
      // MF   
      return intParam[0]*apsi*apsi;
    } else if (intChoice == 2) {
      // 1D-LHY
      return intParam[0]*apsi*apsi+0.5*intParam[1]*apsi;
    } else if (intChoice == 3) {  
      // 2D-LHY   
      return apsi*apsi*(log(apsi*apsi)+1.0);
    } else {  
      printf("Invalid Uintchoice.\n");
      perror(0);
    }

    return 0.0;
} 

double d_Bint2C(double apsi1, double apsi2, const double *intParam,
		const int intChoice, int comp) {
  double n1 = apsi1*apsi1;
  double n2 = apsi2*apsi2;
  if(!intChoice) {
    return 0.0;
  } else if (intChoice == 1 && comp == 1) {
    // MF 
    return intParam[1]*n1;
  } else if (intChoice == 1 && comp == 2) {
    // MF
    return intParam[1]*n2;
    
  } else if (intChoice == 3 && comp == 1) {
    return intParam[1]*n1/intParam[0]-intParam[1]*sqrt(intParam[1])/
      (M_PI*intParam[0]*intParam[0])*intParam[3]*
      n1/(2.0*sqrt(n1/intParam[0]+n2*intParam[0]));
  } else if (intChoice == 3 && comp == 2) {
    return intParam[1]*n2*intParam[0]-intParam[1]*sqrt(intParam[1])*
      intParam[0]*intParam[0]/M_PI*intParam[3]*
      n2/(2.0*sqrt(n1/intParam[0]+n2*intParam[0]));
    
  } else if (intChoice == 4 && comp == 1) {
    return (intParam[0]-intParam[0]*sqrt(intParam[0])/
	    (2.0*M_PI*sqrt(n1+n2)))*n1;
  } else if (intChoice == 4 && comp == 2) {
    return (intParam[0]-intParam[0]*sqrt(intParam[0])/
	    (2.0*M_PI*sqrt(n1+n2)))*n2;
  }
  return 0.0;
}
/*
double complex z_Bint(double complex psi, const double *intParam,
		      const int intChoice) {
  if(!intChoice) {
    return 0.0;
  } else if (intChoice == 1) {
    // MF
    return intParam[0]*psi*psi;
  } else if (intChoice == 2) {
    // 2D-LHY
    return psi*psi*log(psi*conj(psi));
  } else {
    printf("Invalid Uintchoice.\n");
    perror(0);
  }

  return 0.0;
}
*/

void setup_C(coo_t *cooC, double *apsi1, double *apsi2, const int *N,
	     const int dim, const double *intParam, const int intChoice) {
  int size = 1;
  double dapsi;
  for(int i = 0; i < dim; i++)
    size *= N[i];
  for(int i = 0; i < size; i++) {
    dapsi = apsi1[i]*apsi2[i];
    coo_insert(cooC, d_Cint(dapsi, apsi1[i], apsi2[i], intParam, intChoice), 
	       0.0, i, i);
  }

}

double d_Cint(double dapsi, double apsi1, double apsi2, 
	      const double *intParam, const int intChoice) {
  // all need to include a factor of 2.0 times the original derivation
  if(!intChoice) 
    return 0.0;
  else if (intChoice == 1) 
    return 2.0*intParam[2]*dapsi;
  else if (intChoice == 3) {
    return 2.0*(intParam[2]*dapsi-intParam[3]*intParam[1]*sqrt(intParam[1])/M_PI*dapsi/
		(2.0*sqrt(apsi1*apsi1/intParam[0]+apsi2*apsi2*intParam[0])));
  } else if (intChoice == 4) {
    return 2.0*(intParam[1]-intParam[0]-intParam[0]*sqrt(intParam[0])/
	    (2.0*M_PI*sqrt(apsi1*apsi1+apsi2*apsi2)))*apsi1*apsi2;
  }

  return 0.0;
}

void setup_KM1C(coo_t *K, coo_t *diff, coo_t *trap, coo_t *uint,   
		coo_t *B1, const char mode, const int nnz,
		const char flag) {  

  coo_t *work1, *work2;
  work1 = coo_malloc(nnz, trap->rows, trap->cols, flag);
  work2 = coo_malloc(nnz, trap->rows, trap->cols, flag);

  if(mode == 'K') { 
    coo_add(trap, uint, work1);
    //coo_add(work1, B1, work2);
    coo_sub(work1, B1, work2);
    coo_add(diff, work2, K);
  } else if (mode == 'M') { 
    coo_add(trap, uint, work1);
    //coo_sub(work1, B1, work2);
    coo_add(work1, B1, work2);
    coo_add(diff, work2, K);
  } 

  coo_free( work1 );
  coo_free( work2 );
  
}

void setup_KM2C(coo_t *KM, coo_t *diff, coo_t* trap, coo_t *Uint1, coo_t *Uint2,
		coo_t *B1, coo_t *B2, coo_t *C, const char mode, const int nnz,
		const char flag) {

  coo_t *work1, *work2, *work3, *work4, *workM1, *workM2;
  work1 = coo_malloc(nnz, trap->rows, trap->cols, flag);
  work2 = coo_malloc(nnz, trap->rows, trap->cols, flag);
  work3 = coo_malloc(nnz, trap->rows, trap->cols, flag);
  work4 = coo_malloc(nnz, trap->rows, trap->cols, flag);
  //  printf("C: %d\t%d\n",nnz,C->nnz);
  workM1 = coo_malloc(nnz+C->nnz, 2*trap->rows, 2*trap->cols, flag);
  workM2 = coo_malloc(nnz+C->nnz, 2*trap->rows, 2*trap->cols, flag);

  // build A1 and A2
  //  printf("build A1 and A2 %d\t%d\t%d\n",diff->nnz,trap->nnz,work1->nnz);
  coo_add(diff, trap, work1);
  coo_add(work1, Uint1, work2);
  coo_add(work1, Uint2, work3);

  if(mode == 'K') {
    coo_reset(work1);
    coo_sub(work2, B1, work1);
    coo_sub(work3, B2, work4);
    // merge the two by shifting indices of B2 by trap->rows
    // and then use add
    for(int i = 0; i < diff->nnz; i++) {
      work4->col[i] += trap->cols;
      work4->row[i] += trap->rows;
    }
    coo_add(work4, work1, KM);
  } else if(mode == 'M') {
    coo_reset(work1);
    coo_add(work2, B1, work1);
    coo_add(work3, B2, work4);
    // merge the four by shifting indices of 2C and A2+B2
    for(int i = 0; i < diff->nnz; i++) {
      work4->col[i] += trap->cols;
      work4->row[i] += trap->rows;
    }
    for(int i = 0; i < C->nnz; i++)
      C->col[i] += trap->cols;
    //    printf("work1, C, workM1 %d\t%d\t%d\n",work1->nnz, C->nnz, workM1->nnz);
    //    printf("this is a suspicious one\n");
    coo_add(work1, C, workM1);
    for(int i = 0; i < C->nnz; i++) {
      C->col[i] -= trap->cols;
      C->row[i] += trap->cols;
    }
    coo_add(work4, C, workM2);
    //    printf("this is a suspicious one\n");
    coo_add(workM1, workM2, KM);
  }

  coo_free( work1 );
  coo_free( work2 );
  coo_free( work3 );
  coo_free( work4 );
  coo_free( workM1 );
  coo_free( workM2 );

}


void init_eig_vecd(const int N, const int nev, double *vec, double low, double high) {
  
  for(int i = 0;i < nev;i++) {
    for(int j = 0;j < N;j++) {
      vec[j+i*N] = real_rand(low, high);
    }
  }
  
}

void init_eig_vecz(const int N, const int nev, double complex *vec,
		   double low, double high) {
  for(int i = 0; i < nev; i++) {
    for(int j = 0; j < N; j++) {
      vec[j+i*N] = CMPLX(real_rand(low,high),0.0);
      //      vec[j+i*N] = CMPLX(real_rand(low,high),real_rand(low,high));
    } 
  }

}

void init_eig_vec(const int N, const int nev, void *vec,
		  double low, double high, const char flag) {

  if(flag == 'D') {
    init_eig_vecd(N, nev, vec, low, high);
  } else if(flag == 'Z') {
    init_eig_vecz(N, nev, vec, low, high);
  }

}
/*
void add_mode2C(const int dim, const int N[dim], double complex *mode,
		const char wfString[], const char saveString[]) {

  double size;
  for(int i = 0; i < dim; i++) 
    size += N[i];

  double complex *wf1, *wf2;
  wf1 = xmalloc(size*sizeof(double complex));
  wf2 = xmalloc(size*sizeof(double complex));

  load_wf_2C(wf1, wf2, size, wfString);

  if(dim == 1) {
    for(int i = 1; i < N[0]-1; i++) {
      wf1[i] = wf1[i] + mode[i];
      wf2[i] = wf2[i] + mode[i+N];
    }
  }

  safe_free( wf1 );
  safe_free( wf2 );
}
*/
