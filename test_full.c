#include "../LREP/src/splrep.h"
#include "../LREP/src/splopb4dcg.h"
#include <rsb.h>

int main(void) {

  const int type = 2;
  const double shift = 0.0;
  const int rsbTune = 0;

  const int dim = 1;

  int N[dim];
  N[0] = 128;

  double ln[dim];
  ln[0] = 2.0*M_PI;

  const int nev = 8;
  const int subspaceSize = 30;
  const int trapChoice = 1;
  const int intChoice = 2;
  const double tol = 1.0e-8;
  const int maxIter = 100;

  rsb_time_t dt;

  double param[5] = { 0.0 };
  double intParam[5] = { 0.0 };

  double partN, g, g12, f1, f2;

  int numsteps = 10;
  double *mu;
  mu = xmalloc(numsteps*sizeof(double));

  char wfFilename[1024];
  char specFilename[1024];
  char vecFilename[1024];
  char enFilename[1024];

  int counter = 0;
  int fileCounter = 0;
  double buffer = 0.0;

  
  // open file for chemical potential
  snprintf(enFilename, sizeof(enFilename), "chemPotvar12_.dat");
  FILE* inStream = fopen(enFilename,"r");
  if(inStream) {
    while(fscanf(inStream, "%lf", &buffer) && (counter < numsteps)) {
      mu[counter] = buffer;
      counter ++;
    }
    fclose(inStream);
  } else {
    // Provides some error diagnostic
    fprintf(stderr, "Could not open file for chemical potential.\n");
    perror(0);
    errno = 0;
  }
  
  snprintf(specFilename, sizeof(specFilename), "spectrum_exact_var12.dat");
  FILE* outStream = fopen(specFilename, "w");
  if(!outStream) {
    printf("Could not open file to save eigenvalues.\n");
    return 0;
  }
  
  printf("also reached this point.\n");
  for (int j = 0; j < numsteps; j++) {
    
    snprintf(vecFilename, sizeof(vecFilename), "vec_exact_var12_%d.dat",j);
    FILE* outStream1 = fopen(vecFilename, "w");
    if(!outStream1) {
      printf("Could not open file to save eigenvectors.\n");
      return 0;
    }
    
    g = 2.5/numsteps*j;
    g12 = -0.9*g;
    partN = 12.0;
    f1 = g+g12;
    f2 = g-g12;
    
    intParam[0] = partN/2.0*f1;
    intParam[1] = -sqrt(partN)/(M_PI*sqrt(2.0)*2.0)*(f1*sqrt(f1)+f2*sqrt(f2));
    
    struct sp_lrep_t *LREP;
    
    snprintf(wfFilename, sizeof(wfFilename), "filepsi_1D_gvar12_%d.dat", j);
    dt = -rsb_time();
    LREP = splrep_setup(wfFilename, N, ln, nev, subspaceSize, dim,
			mu[j], trapChoice, param, intChoice, intParam, type, shift, maxIter,
			tol, rsbTune, DATA_TYPE_COMPLEX);
    dt += rsb_time();
    printf("Setting up LREP took %lfs.\n",dt);
    
    dt = -rsb_time();
    sp_lopb4dcg(LREP);
    dt += rsb_time();
    printf("LOPB4DCG took %lfs.\n",dt);
    
    for(int i = 0; i < nev; i++) {
      fprintf(outStream, "%.7e\t%.7e\t",
	      creal(LREP->eValSort[i]),cimag(LREP->eValSort[i]));
    }
    fprintf(outStream,"\n");

    for(int i = 0; i < LREP->nev; i++) {
      for(int k = 0; k < LREP->size; k++) {
	fprintf(outStream1, "%.7e\t%.7e\t%.7e\t%.7e\n", creal(LREP->X[k*subspaceSize+i]),
		cimag(LREP->X[k*subspaceSize+i]),creal(LREP->Y[k*subspaceSize+i]),
		cimag(LREP->Y[k*subspaceSize+i]));
      }
    }
    
    splrep_free( LREP );
    fclose(outStream1);

    /*    for(int i = 0; i < LREP->nev; i++) {                                                                           
      for(int k = 0; j < LREP->size; k++) {           
	fprintf(outStream1, "%.7e\t%.7e\t%.7e\t%.7e\n", creal(LREP->X[k*subspaceSize+i]),                          
		cimag(LREP->X[k*subspaceSize+i]));//,creal(LREP->Y[k*subspaceSize+i]),
		//cimag(LREP->Y[k*subspaceSize+i]));                                                                 
      }                                                                                                            
    }   
    */
  }
  
  fclose(outStream);
  
  
  safe_free( mu );
  

  return 0;
}
