
#include "splopb4dcg.h"
#include <rsb.h>

void sp_get_search_dird(struct sp_lrep_t *LREP) {

  double num1, num2, denom;
  const int n = LREP->size;
  const int sizeSub = LREP->sizeSub;
  double *work;
  work = xmalloc(LREP->size*sizeof(double));

  double *transX, *transY;
  transX = xmalloc(n*sizeSub*sizeof(double));
  transY = xmalloc(n*sizeSub*sizeof(double));
  matrix_trans(n, sizeSub, (double *)LREP->X, transX, LREP->flag);
  matrix_trans(n, sizeSub, (double *)LREP->Y, transY, LREP->flag);

  // calculate the first rho as initial start of the calculation
  // calculate M*X from x^T*M*x, then the dot product

  for(int i = 0; i < sizeSub; i++) {
    rsb_SPMV(LREP->K, &transX[i*n], work, LREP->flag);
    num1 = ddot_product(n, &transX[i*n], work);

    rsb_SPMV(LREP->M, &transY[i*n], work, LREP->flag);
    num2 = ddot_product(n, &transY[i*n], work);

    // calculate the denominator
    denom = ddot_product(n, &transX[i*n], &transY[i*n]);
    denom = 2*fabs(denom);
    // this will throw some error
    ((double *)LREP->eigVal)[i+i*sizeSub] = (num1 + num2)/denom;
  }

  safe_free( work );
  safe_free( transX );
  safe_free( transY );

  
}

void sp_get_search_dirz(struct sp_lrep_t *LREP) {

  double complex num1, num2, denom;
  const int n = LREP->size;
  const int sizeSub = LREP->sizeSub;
  double complex *work;
  work = xmalloc(LREP->size*sizeof(double complex));

  double complex *transX, *transY;
  transX = xmalloc(n*sizeSub*sizeof(double complex));
  transY = xmalloc(n*sizeSub*sizeof(double complex));
  matrix_trans(n, sizeSub, (double complex *)LREP->X, transX, LREP->flag);
  matrix_trans(n, sizeSub, (double complex *)LREP->Y, transY, LREP->flag);

  // calculate the first rho as initial start of the calculation
  // calculate M*X from x^T*M*x, then the dot product

  for(int i = 0; i < sizeSub; i++) {
    rsb_SPMV(LREP->K, &transX[i*n], work, LREP->flag);
    num1 = zdot_product(n, &transX[i*n], work);

    rsb_SPMV(LREP->M, &transY[i*n], work, LREP->flag);
    num2 = zdot_product(n, &transX[i*n], work);

    // calculate the denominator
    denom = zdot_product(n, &transX[i*n], &transX[i*n]);
    denom = 2*cabs(denom);
    ((double complex *)LREP->eigVal)[i+i*sizeSub] = (num1 + num2)/denom;
  }

  safe_free( work );
  safe_free( transX );
  safe_free( transY );

  
}

void sp_get_search_dir(struct sp_lrep_t *LREP) {
  switch(LREP->flag) {
  case 'D':
    sp_get_search_dird(LREP);
    break;

  case 'Z':
    sp_get_search_dirz(LREP);
    break;
  }
}

void sp_setup_P(struct sp_lrep_t *LREP) {    
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;   
  
  // calculate Y*diag(search_dir)
  gemm_NN(n, sizeSub, sizeSub, LREP->Y, LREP->eigVal, LREP->P, LREP->flag);
  
  // calculate K*X-Y*diag(search_dir) 
  rsb_SPMM_sub(LREP->K, LREP->X, LREP->P, sizeSub, LREP->flag);      
}  

void sp_setup_Q(struct sp_lrep_t *LREP) {    
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;   
   
  // calculate X*diag(search_dir)
  gemm_NN(n, sizeSub, sizeSub, LREP->X, LREP->eigVal, LREP->Q, LREP->flag);
  
  // calculate M*Y-X*diag(search_dir) 
  rsb_SPMM_sub(LREP->M, LREP->Y, LREP->Q, sizeSub, LREP->flag);   
} 

void sp_apply_precondd(struct sp_lrep_t *LREP) { 
   
  int type = LREP->spPrecond->type;  
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;
  double *work1 = NULL, *work2 = NULL;
  work1 = xcalloc(n*sizeSub, sizeof(double));
  work2 = xcalloc(n*sizeSub, sizeof(double));

  switch(type) {
  case 0:
    // do nothing as not preconditioner chosen
    break;

  case 1:
    // diagonal preconditioner
    // multiply Mdiag and Kdiag on P and Q respectively   
    rsb_SPMM(LREP->spPrecond->KDiag, LREP->P, work1, sizeSub, LREP->flag); 
    rsb_SPMM(LREP->spPrecond->MDiag, LREP->Q, work2, sizeSub, LREP->flag); 
   
    copy_vec(n*sizeSub, work1, LREP->P, LREP->flag);    
    copy_vec(n*sizeSub, work2, LREP->Q, LREP->flag);    
  
    break;

  case 2:
    // ILU(0)
    // first apply Preconditioner on P, which then turns to Q  
    rsb_SPSM(LREP->spPrecond->MPrecond[0], LREP->P, work1, sizeSub, LREP->flag);
    rsb_SPSM(LREP->spPrecond->KPrecond[0], LREP->Q, work2, sizeSub, LREP->flag);
   
    rsb_SPSM(LREP->spPrecond->MPrecond[1], work1, LREP->Q, sizeSub, LREP->flag);
    rsb_SPSM(LREP->spPrecond->KPrecond[1], work2, LREP->P, sizeSub, LREP->flag);

    break;

  case 3:
    // block_conjugate gradient
    sp_block_conj_gradd(LREP->K, (double *)LREP->P, work1, n, sizeSub,
		       LREP->spPrecond->KDiag, 20, 1.0e-4);
    sp_block_conj_gradd(LREP->M, (double *)LREP->Q, work2, n, sizeSub,
		       LREP->spPrecond->MDiag, 20, 1.0e-4);

    copy_vec(n*sizeSub, work1, LREP->P, LREP->flag);
    copy_vec(n*sizeSub, work2, LREP->Q, LREP->flag);

    break;

  default:
    fprintf(stderr, "Preconditioner type not recognized.\n");
    break;
  }

  safe_free( work1 );
  safe_free( work2 );
      
}

void sp_apply_precondz(struct sp_lrep_t *LREP) { 
   
  int type = LREP->spPrecond->type;  
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;
  double complex *work1 = NULL, *work2 = NULL;
  work1 = xcalloc(n*sizeSub, sizeof(double complex));
  work2 = xcalloc(n*sizeSub, sizeof(double complex));

  switch(type) {
  case 0:
    // do nothing as not preconditioner chosen
    break;

  case 1:
    // diagonal preconditioner
    // multiply Mdiag and Kdiag on P and Q respectively   
    rsb_SPMM(LREP->spPrecond->KDiag, LREP->P, work1, sizeSub, LREP->flag); 
    rsb_SPMM(LREP->spPrecond->MDiag, LREP->Q, work2, sizeSub, LREP->flag); 
   
    copy_vec(n*sizeSub, work1, LREP->P, LREP->flag);    
    copy_vec(n*sizeSub, work2, LREP->Q, LREP->flag);    
  
    break;

  case 2:
    // ILU(0)
    // first apply Preconditioner on P, which then turns to Q  
    rsb_SPSM(LREP->spPrecond->MPrecond[0], LREP->P, work1, sizeSub, LREP->flag);
    rsb_SPSM(LREP->spPrecond->KPrecond[0], LREP->Q, work2, sizeSub, LREP->flag);
   
    rsb_SPSM(LREP->spPrecond->MPrecond[1], work1, LREP->Q, sizeSub, LREP->flag);
    rsb_SPSM(LREP->spPrecond->KPrecond[1], work2, LREP->P, sizeSub, LREP->flag);
   
    break;

  case 3:
    // block_conjugate gradient
    sp_block_conj_gradz(LREP->K, LREP->P, work1, n, sizeSub,
		       LREP->spPrecond->KDiag, 20, 1.0e-4);
    sp_block_conj_gradz(LREP->M, LREP->Q, work2, n, sizeSub,
		       LREP->spPrecond->MDiag, 20, 1.0e-4);

    copy_vec(n*sizeSub, work1, LREP->P, LREP->flag);
    copy_vec(n*sizeSub, work2, LREP->Q, LREP->flag);
   
    break;

  default:
    fprintf(stderr, "Preconditioner type not recognized.\n");
    break;
  }

  safe_free( work1 );
  safe_free( work2 );
      
}

void sp_apply_precond(struct sp_lrep_t *LREP) {

  switch(LREP->flag) {
  case 'D':
    sp_apply_precondd(LREP);
    break;

  case 'Z':
    sp_apply_precondz(LREP);
    break;
  }
}


void sp_setup_U(struct sp_lrep_t *LREP) {
  // could also just do the whole thing with copy_vec??
  // can also set the initial pointers to U and V
  // to avoid repeated accesses
  const int n = LREP->size;
  const int sizeSub = LREP->sizeSub;
  const int iter = LREP->iter;

  switch(iter) {
    
  case 0:

#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
      /*for(int j = 0; j < sizeSub*mult; j++) {
	
	if(j < sizeSub) {
	  LREP->U[j + i*sizeSub*mult] = LREP->X[j+i*sizeSub];
	} else {
	  LREP->U[j + i*sizeSub*mult] = LREP->P[j-sizeSub+i*sizeSub];
	}
	}*/
      for(int j = 0; j < sizeSub; j++) {
	LREP->U[j + i*sizeSub] = LREP->X[j + i*sizeSub];
	LREP->U[j+sizeSub + i*sizeSub] = LREP->P[j+i*sizeSub];
      }
    }
    break;

  default:

#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < sizeSub; j++) {
	LREP->U[j + i*sizeSub] = LREP->X[j+i*sizeSub];
	LREP->U[j+sizeSub + i*sizeSub] = LREP->X1[j+i*sizeSub];
	LREP->U[j+2*sizeSub + i*sizeSub] = LREP->P[j+i*sizeSub];
      }
    }
    break;
  }
}
  /*
      for(int j = 0; j < sizeSub*mult; j++) {
	if(j < sizeSub) {
	  LREP->U[j+i*sizeSub*mult] = LREP->X[j+i*sizeSub];
	} else if (j < 2*sizeSub && j >= sizeSub) {
	  LREP->U[j+i*sizeSub*mult] = LREP->X1[j-sizeSub+i*sizeSub];
	} else {
	  LREP->U[j+i*sizeSub*mult] = LREP->P[j-2*sizeSub+i*sizeSub];
	}
      }
    }
  */
/*  
  for(int i = 0; i < n; i++) {   
    for(int j = 0; j < sizeSub*mult; j++) { 
      if(LREP-> iter == 0) {
	if(j < sizeSub) {    
	  LREP->U[j+i*sizeSub*mult] = LREP->X[j+i*sizeSub]; 
	} else {  
	  LREP->U[j+i*sizeSub*mult] = LREP->P[j-sizeSub+i*sizeSub];    
	}    
      } else {    
	if(j < sizeSub) {    
	  LREP->U[j+i*sizeSub*mult] = LREP->X[j+i*sizeSub]; 
	} else if (j < 2*sizeSub && j >= sizeSub) {    
	  LREP->U[j+i*sizeSub*mult] = LREP->X1[j-sizeSub+i*sizeSub];   
	} else {  
	  LREP->U[j+i*sizeSub*mult] = LREP->P[j-2*sizeSub+i*sizeSub];  
	}    
      } 
    }   
  }   
}
*/
void sp_setup_V(struct sp_lrep_t *LREP) {
  // could also just replace with copy_vec?
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;
  const int iter = LREP->iter;

  switch(iter) {
  case 0:
#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < sizeSub; j++) {
	LREP->V[j + i*sizeSub] = LREP->Y[j+i*sizeSub];
	LREP->V[j+sizeSub + i*sizeSub] = LREP->Q[j+i*sizeSub];
      }
    }

    break;

  default:
#pragma omp parallel for collapse(2)
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < sizeSub; j++) {
	LREP->V[j + i*sizeSub] = LREP->Y[j+i*sizeSub];
	LREP->V[j+sizeSub + i*sizeSub] = LREP->Y1[j+i*sizeSub];
	LREP->V[j+2*sizeSub + i*sizeSub] = LREP->Q[j+i*sizeSub];
      }
    }

    break;
  }
}

  /*
  for(int i = 0; i < n; i++) {   
    for(int j = 0; j < sizeSub*mult; j++) { 
      if(LREP-> iter == 0) {
	if(j < sizeSub) {    
	  LREP->V[j+i*sizeSub*mult] = LREP->Y[j+i*sizeSub]; 
	} else {  
	  LREP->V[j+i*sizeSub*mult] = LREP->Q[j-sizeSub+i*sizeSub];    
	}    
      } else {    
	if(j < sizeSub) {    
	  LREP->V[j+i*sizeSub*mult] = LREP->Y[j+i*sizeSub]; 
	} else if (j < 2*sizeSub && j >= sizeSub) {    
	  LREP->V[j+i*sizeSub*mult] = LREP->Y1[j-sizeSub+i*sizeSub];   
	} else {  
	  LREP->V[j+i*sizeSub*mult] = LREP->Q[j-2*sizeSub+i*sizeSub];  
	}    
      } 
    }   
    }
    }*/

void sp_setup_W(struct sp_lrep_t *LREP) {    
  const int n = LREP->size;
  const int mult = (0 == LREP->iter) ? 2 : 3;
  const int sizeSub = LREP->sizeSub*mult;   
  
  gemm_TN(sizeSub, sizeSub, n, LREP->U, LREP->V, LREP->W, LREP->flag);   
}

void sp_setup_Hsrd(struct sp_lrep_t *LREP) {  
  // I am not so super sure about the sizes for rsb_SPMM here  
  const int n = LREP->size;
  const int mult = (0 == LREP->iter) ? 2 : 3;
  const int sizeSub = LREP->sizeSub*mult;   
  
  double *work, *work1, *work2, *work3;    
  work = xmalloc(n*sizeSub*sizeof(double)); 
  work1 = xmalloc(sizeSub*sizeSub*sizeof(double));
  work2 = xmalloc(n*sizeSub*sizeof(double));
  work3 = xmalloc(sizeSub*sizeSub*sizeof(double));
  
  // Set all entries in Hsr to 0 to get rid of  
  // diagonalization from previous iteration    
  memset(LREP->Hsr, 0, 4*sizeSub*sizeSub*sizeof(double));   
  
  // Compute top right block matrix first
  rsb_SPMM(LREP->K, LREP->U, work, sizeSub, LREP->flag);
  gemm_TN(sizeSub, sizeSub, n, LREP->U, work, work1, LREP->flag);

#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {   
    for(int j = 0; j < sizeSub; j++) { 
      LREP->Hsr[j+sizeSub*(2*i+1)] = work1[j+i*sizeSub];    
    }   
  }

  // Compute bottom left block matrix 
  // Assumes that W is inv(W)    
  gemm_NN(n, sizeSub, sizeSub, LREP->V, LREP->W, work, LREP->flag);   
  rsb_SPMM(LREP->M, work, work2, sizeSub, LREP->flag);  

  gemm_TN(sizeSub, sizeSub, n, LREP->V, work2, work3, LREP->flag);    
  gemm_TN(sizeSub, sizeSub, sizeSub, LREP->W, work3, work1, LREP->flag);

#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {   
    for(int j = 0; j < sizeSub; j++) {  
      LREP->Hsr[2*sizeSub*sizeSub+j+2*sizeSub*i] = work1[j+i*sizeSub];  
    }   
  }

  safe_free( work );   
  safe_free( work1 );  
  safe_free( work2 );  
  safe_free( work3 ); 
  
  
}

void sp_setup_Hsrz(struct sp_lrep_t *LREP) {
  const int n = LREP->size;
  const int mult = (0 == LREP->iter) ? 2 : 3;
  const int sizeSub = LREP->sizeSub*mult;

  double complex *work, *work1, *work2, *work3;
  work = xmalloc(n*sizeSub*sizeof(double complex));
  work1 = xmalloc(sizeSub*sizeSub*sizeof(double complex));
  work2 = xmalloc(n*sizeSub*sizeof(double complex));
  work3 = xmalloc(sizeSub*sizeSub*sizeof(double complex));

  // Set all entries in Hsr to 0 to get rid of
  // diagonalization from previous iteration
  memset(LREP->Hsr, 0, 4*sizeSub*sizeSub*sizeof(double complex));

  // Compute top right block matrix first
  rsb_SPMM(LREP->K, LREP->U, work, sizeSub, LREP->flag);
  gemm_TN(sizeSub, sizeSub, n, LREP->U, work, work1, LREP->flag);
#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {
    for(int j = 0; j < sizeSub; j++) {
      LREP->Hsr[j+sizeSub*(2*i+1)] = work1[j+i*sizeSub];
    }
  }

  // Compute bottom left block matrix
  // Assumes that W is inv(W)
  gemm_NN(n, sizeSub, sizeSub, LREP->V, LREP->W, work, LREP->flag);
  rsb_SPMM(LREP->M, work, work2, sizeSub, LREP->flag);
  gemm_TN(sizeSub, sizeSub, n, LREP->V, work2, work3, LREP->flag);
  gemm_TN(sizeSub, sizeSub, sizeSub, LREP->W, work3, work1, LREP->flag);
#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {
    for(int j = 0; j < sizeSub; j++) {
      LREP->Hsr[2*sizeSub*sizeSub+j+2*sizeSub*i] = work1[j+i*sizeSub];
    }
  }

  safe_free( work );
  safe_free( work1 );
  safe_free( work2 );
  safe_free( work3 );
}

void sp_setup_Hsr(struct sp_lrep_t *LREP) {

  switch(LREP->flag) {
  case 'Z':
    sp_setup_Hsrz(LREP);
    break;

  default:
    sp_setup_Hsrd(LREP);
    break;

  }
}
 

void sp_sort_eigd(struct sp_lrep_t *LREP) {   
  const int sizeSub = LREP->sizeSub;
  const int mult = (0 == LREP->iter) ? 4 : 6;
  
  int arrSorted[sizeSub];    
  double curr = 0.0;   
  double last_min = 0.0;    
  int idx = 0;    
  double tol = 1e-6;   
  
  double eVecT[mult*sizeSub*mult*sizeSub];   
  double eVecTSort[sizeSub*sizeSub*mult];    
  matrix_trans(mult*sizeSub, mult*sizeSub, LREP->eVecsr, eVecT, LREP->flag);  
  
  for(int i = 0; i < sizeSub; i++) {   
    curr = INFINITY;   
    for(int j = 0; j < mult*sizeSub; j++) { 
      if(creal(LREP->eValsr[j]) < tol ||
	 creal(LREP->eValsr[j]) <= last_min ||    
	 creal(LREP->eValsr[j])/last_min - 1.0 < tol) {
	continue; 
      } else {    
	if(creal(LREP->eValsr[j]) < curr) {  
	  curr = LREP->eValsr[j];
	  idx = j;
	} else {  
	  continue;    
	}    
      } 
    }   
    arrSorted[i] = idx;
    LREP->eValSort[i] = curr;    
    last_min = curr;   
  }

#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {   
    for(int j = 0; j < mult*sizeSub; j++) { 
      eVecTSort[j+i*mult*sizeSub] = eVecT[j+arrSorted[i]*mult*sizeSub];    
    }   
  }
  matrix_trans(sizeSub, mult*sizeSub, eVecTSort, LREP->eVecSort, LREP->flag); 
      
}

void sp_sort_eigz(struct sp_lrep_t *LREP) {

  // for complex case this works a bit differently
  // first have to sort all eigenvalues (or 2*sizeSub)
  // according to amplitude
  // then go through those pairs and pick the positive one
  // (and somehow keep the index of everything for eigenvectors)
  
  const int sizeSub = LREP->sizeSub;
  const int mult = (0 == LREP->iter) ? 4 : 6;
  int arrSort[sizeSub];
  int arrSortAmp[mult*sizeSub];

  double complex *eVecT, *eVecTSort;
  eVecT = xmalloc(mult*sizeSub*mult*sizeSub*sizeof(double complex));
  eVecTSort = xmalloc(mult*sizeSub*sizeSub*sizeof(double complex));
  matrix_trans(mult*sizeSub, mult*sizeSub, LREP->eVecsr, eVecT, LREP->flag);

  double tol = 1.0e-13;
  double complex tempz = -1.0;
  int tempi = -1;
  int swapped = 0;

  // bubble sort eigenvalues according to their complex amplitude

#pragma omp parallel for
  for(int i = 0; i < mult*sizeSub; i++) {
    arrSortAmp[i] = i;
  }

  //#pragma omp parallel for collapse(2)
  for(int i = 0; i < mult*sizeSub-1; i++) {
    for(int j = 0; j < mult*sizeSub-i-1; j++) {
      if(cabs(LREP->eValsr[j]) > cabs(LREP->eValsr[j+1])) {
	tempz = LREP->eValsr[j];
	LREP->eValsr[j] = LREP->eValsr[j+1];
	LREP->eValsr[j+1] = tempz;

	tempi = arrSortAmp[j];
	arrSortAmp[j] = arrSortAmp[j+1];
	arrSortAmp[j+1] = tempi;

	swapped = 1;
      }
    }

    if(!swapped)
      break;
  }

  // now that I have sorted it according to amplitude, need to go through it
  // and check which one of them are the positive ones
  // three cases: creal and cimag positive
  // creal pos and cimag neg
  // creal neg and cimag pos

  int counter = 0;
  int prevCounter = 0;
  int sel = 0;
  //  for(int i = 0; i < mult*sizeSub; i++)
  //    printf("%.5e %.5e\t",creal(LREP->eValsr[i]),cimag(LREP->eValsr[i]));
  //  printf("\n");
  for(int i = 0; i < sizeSub; i++) {
    prevCounter = counter;
    counter = 0;
    sel = 0;
    if(creal(LREP->eValsr[2*i]) > tol ||
       cimag(LREP->eValsr[2*i]) > tol) {
      counter++;
      sel = 1;
    }
    if(creal(LREP->eValsr[2*i+1]) > tol ||
       cimag(LREP->eValsr[2*i+1]) > tol) {
      counter++;
      sel = 2;
    }

    //    printf("%d\t", counter);

    if(counter == 1 && sel == 1) {
      arrSort[i] = arrSortAmp[2*i];
      LREP->eValSort[i] = LREP->eValsr[2*i];
    }
    if(counter == 1 && sel == 2) {
      arrSort[i] = arrSortAmp[2*i+1];
      LREP->eValSort[i] = LREP->eValsr[2*i+1];
    }

    if(counter == 2 && !prevCounter && (i > 0)) {
      arrSort[i] = arrSortAmp[2*i];
      LREP->eValSort[i] = LREP->eValsr[2*i];
      arrSort[i-1] = arrSortAmp[2*i+1];
      LREP->eValSort[i-1] = LREP->eValsr[2*i+1];
    }
    if(!counter && prevCounter == 2) {
      arrSort[i] = arrSortAmp[2*i-1];
      LREP->eValSort[i] = LREP->eValsr[2*i-1];
      arrSort[i-1] = arrSortAmp[2*i-2];
      LREP->eValSort[i-1] = LREP->eValsr[2*i-2];
    }

    if((counter == 2 && prevCounter) || (counter == 2 && !i)) {
      if(fabs(creal(LREP->eValsr[2*i])) > fabs(cimag(LREP->eValsr[2*i]))) {
	if(creal(LREP->eValsr[2*i]) > creal(LREP->eValsr[2*i+1])) {
	  arrSort[i] = arrSortAmp[2*i];
	  LREP->eValSort[i] = LREP->eValsr[2*i];
	} else {
	  arrSort[i] = arrSortAmp[2*i+1];
	  LREP->eValSort[i] = LREP->eValsr[2*i+1];
	}
      } else {
	if(cimag(LREP->eValsr[2*i]) > creal(LREP->eValsr[2*i+1])) {
	  arrSort[i] = arrSortAmp[2*i];
	  LREP->eValSort[i] = LREP->eValsr[2*i];
	} else {
	  arrSort[i] = arrSortAmp[2*i+1];
	  LREP->eValSort[i] = LREP->eValsr[2*i+1];
	}
      }
    }
    //    printf("%d\t",counter);
  }
  //  printf("\n");
  //  print_array(sizeSub,1,arrSort);

#pragma omp parallel for collapse(2)
  for(int i = 0; i < sizeSub; i++) {
    for(int j = 0; j < mult*sizeSub; j++) {
      eVecTSort[j+i*mult*sizeSub] = eVecT[j+arrSort[i]*mult*sizeSub];
    }
  }

  matrix_trans(sizeSub, mult*sizeSub, eVecTSort, LREP->eVecSort, LREP->flag);

  safe_free( eVecT );
  safe_free( eVecTSort );

  
}

void sp_sort_eig(struct sp_lrep_t *LREP) {
  switch(LREP->flag) {
  case 'Z':
    sp_sort_eigz(LREP);
    break;

  default:
    sp_sort_eigd(LREP);
    break;
  }
}

void sp_split_eig_vec(struct sp_lrep_t *LREP) {    
  const int sizeSub = LREP->sizeSub;
  const int mult = (0 == LREP->iter) ? 4 : 6;

#pragma omp parallel for
  for(int i = 0; i < sizeSub*sizeSub*mult/2; i++) {
    LREP->Ysr[i] = LREP->eVecSort[i];
    LREP->Xsr[i] = LREP->eVecSort[i + sizeSub*sizeSub*mult/2];
  }
}
    /*
    for(int i = 0; i < sizeSub*sizeSub*mult; i++) {
    if(i < sizeSub*sizeSub*mult/2) {    
      LREP->Ysr[i] = LREP->eVecSort[i];    
    } else { 
      LREP->Xsr[i-sizeSub*sizeSub*mult/2] = LREP->eVecSort[i];   
    }   
    }
}*/

void sp_compute_eig_vec(struct sp_lrep_t *LREP) {  
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;
  const int mult = (0 == LREP->iter) ? 2 : 3;
   
  double *workd = NULL;   
  double complex *workz = NULL;

  switch(LREP->flag) {
  case 'Z':
    workz = xmalloc(mult*sizeSub*sizeSub*sizeof(double complex));
    // first do X
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->U, LREP->Xsr, LREP->X1, LREP->flag);
    // then Y
    gemm_NN(mult*sizeSub, sizeSub, mult*sizeSub, LREP->W, LREP->Ysr, workz, LREP->flag);
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->V, workz, LREP->Y1, LREP->flag);
    
    safe_free( workz );

    break;

  default:
    workd = xmalloc(mult*sizeSub*sizeSub*sizeof(double));
    // first do X   
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->U, LREP->Xsr, LREP->X1, LREP->flag);  
    // then Y  
    gemm_NN(mult*sizeSub, sizeSub, mult*sizeSub, LREP->W, LREP->Ysr, workd, LREP->flag); 
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->V, workd, LREP->Y1, LREP->flag);
    
    safe_free( workd );
    break;
  }
}
/*


  if(LREP->flag == 'D') {
    workd = xmalloc(mult*sizeSub*sizeSub*sizeof(double));
    // first do X   
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->U, LREP->Xsr, LREP->X1, LREP->flag);  
    // then Y  
    gemm_NN(mult*sizeSub, sizeSub, mult*sizeSub, LREP->W, LREP->Ysr, workd, LREP->flag); 
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->V, workd, LREP->Y1, LREP->flag);
    
    safe_free( workd );
    
  } else if(LREP->flag == 'Z') {
    workz = xmalloc(mult*sizeSub*sizeSub*sizeof(double complex));
    // first do X
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->U, LREP->Xsr, LREP->X1, LREP->flag);
    // then Y
    gemm_NN(mult*sizeSub, sizeSub, mult*sizeSub, LREP->W, LREP->Ysr, workz, LREP->flag);
    gemm_NN(n, sizeSub, mult*sizeSub, LREP->V, workz, LREP->Y1, LREP->flag);
    
    safe_free( workz );
  }
  }*/


void sp_normalize_eig_vec(struct sp_lrep_t *LREP) {
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;

  double norm, normx, normy;
   
  double *transXd = NULL, *transYd = NULL;  
  double complex *transXz = NULL, *transYz = NULL;

  switch(LREP->flag) {
  case 'Z':
    transXz = xmalloc(n*sizeSub*sizeof(double complex));
    transYz = xmalloc(n*sizeSub*sizeof(double complex));
    
    matrix_trans(n, sizeSub, LREP->X1, transXz, LREP->flag);
    matrix_trans(n, sizeSub, LREP->Y1, transYz, LREP->flag);

    for(int i = 0; i < sizeSub; i++) {
      normx = asum(n, &transXz[i*n], LREP->flag);
      normy = asum(n, &transYz[i*n], LREP->flag);

      norm = normx + normy;

      scale_vec(n, &transXz[i*n], 1.0/norm, LREP->flag);
      scale_vec(n, &transYz[i*n], 1.0/norm, LREP->flag);
    }

    matrix_trans(sizeSub, n, transXz, LREP->X1, LREP->flag);
    matrix_trans(sizeSub, n, transYz, LREP->Y1, LREP->flag);

    safe_free( transXz );
    safe_free( transYz );
    break;

  default:
    transXd = xmalloc(n*sizeSub*sizeof(double));
    transYd = xmalloc(n*sizeSub*sizeof(double));
    
    matrix_trans(n, sizeSub, LREP->X1, transXd, LREP->flag);
    matrix_trans(n, sizeSub, LREP->Y1, transYd, LREP->flag);

    for(int i = 0; i < sizeSub; i++) {   
      normx = asum(n, &transXd[i*n], LREP->flag);
      normy = asum(n, &transYd[i*n], LREP->flag);
      
      norm = normx + normy;   
      
      scale_vec(n, &transXd[i*n], 1.0/norm, LREP->flag);    
      scale_vec(n, &transYd[i*n], 1.0/norm, LREP->flag);    
    }
   
    matrix_trans(sizeSub, n, transXd, LREP->X1, LREP->flag);
    matrix_trans(sizeSub, n, transYd, LREP->Y1, LREP->flag);
    
    safe_free( transXd ); 
    safe_free( transYd );
    break;
  }
}

/*   

  if(LREP->flag == 'D') {

    transXd = xmalloc(n*sizeSub*sizeof(double));
    transYd = xmalloc(n*sizeSub*sizeof(double));
    
    matrix_trans(n, sizeSub, LREP->X1, transXd, LREP->flag);
    matrix_trans(n, sizeSub, LREP->Y1, transYd, LREP->flag);

    for(int i = 0; i < sizeSub; i++) {   
      normx = asum(n, &transXd[i*n], LREP->flag);
      normy = asum(n, &transYd[i*n], LREP->flag);
      
      norm = normx + normy;   
      
      scale_vec(n, &transXd[i*n], 1.0/norm, LREP->flag);    
      scale_vec(n, &transYd[i*n], 1.0/norm, LREP->flag);    
    }
   
    matrix_trans(sizeSub, n, transXd, LREP->X1, LREP->flag);
    matrix_trans(sizeSub, n, transYd, LREP->Y1, LREP->flag);
    
    safe_free( transXd ); 
    safe_free( transYd ); 

  } else if(LREP->flag == 'Z') {

    transXz = xmalloc(n*sizeSub*sizeof(double complex));
    transYz = xmalloc(n*sizeSub*sizeof(double complex));
    
    matrix_trans(n, sizeSub, LREP->X1, transXz, LREP->flag);
    matrix_trans(n, sizeSub, LREP->Y1, transYz, LREP->flag);

    for(int i = 0; i < sizeSub; i++) {
      normx = asum(n, &transXz[i*n], LREP->flag);
      normy = asum(n, &transYz[i*n], LREP->flag);

      norm = normx + normy;

      scale_vec(n, &transXz[i*n], 1.0/norm, LREP->flag);
      scale_vec(n, &transYz[i*n], 1.0/norm, LREP->flag);
    }

    matrix_trans(sizeSub, n, transXz, LREP->X1, LREP->flag);
    matrix_trans(sizeSub, n, transYz, LREP->Y1, LREP->flag);

    safe_free( transXz );
    safe_free( transYz );
  }
      
  }*/

void sp_get_residual_normd(struct sp_lrep_t *LREP) {
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;   
   
  double nominator, denominator; 
  double *work1, *work2;    
  work1 = xmalloc(n*sizeof(double));  
  work2 = xmalloc(n*sizeof(double));  
   
  double *transX, *transY;  
  transX = xmalloc(n*sizeSub*sizeof(double));    
  transY = xmalloc(n*sizeSub*sizeof(double));    
   
  matrix_trans(n, sizeSub, LREP->X1, transX, LREP->flag);
  matrix_trans(n, sizeSub, LREP->Y1, transY, LREP->flag);
   
  // caluclate denominator first 
  // one norm of K&M already calculated upon setup   
  for(int i = 0; i < sizeSub; i++) {   
    denominator = LREP->oneNorm + LREP->eValSort[i]; 
   
    // calculate nominator  
    rsb_SPMV(LREP->K, transX+i*n, work1, LREP->flag);  

    copy_vec(n, transY+i*n, work2, LREP->flag);   
    scale_vec(n, work2, LREP->eValSort[i], LREP->flag);
    sub_vec(n, work1, work2, LREP->flag);    
    nominator = asum(n, work1, LREP->flag);
   
    rsb_SPMV(LREP->M, transY+i*n, work1, LREP->flag);

    copy_vec(n, transX+i*n, work2, LREP->flag);   
    scale_vec(n, work2, LREP->eValSort[i], LREP->flag);
    sub_vec(n, work1, work2, LREP->flag);    
    nominator += asum(n, work1, LREP->flag);
   
    LREP -> resNorm[i] = nominator/denominator; 
  }
   
  safe_free( work1 );  
  safe_free( work2 );  
  safe_free( transX );
  safe_free( transY );
      
}

void sp_get_residual_normz(struct sp_lrep_t *LREP) {
  const int n = LREP->size;
  const int sizeSub = LREP->sizeSub;

  double complex nominator, denominator;
  double complex *work1, *work2;
  work1 = xmalloc(n*sizeof(double complex));
  work2 = xmalloc(n*sizeof(double complex));

  double complex *transX, *transY;
  transX = xmalloc(n*sizeSub*sizeof(double complex));
  transY = xmalloc(n*sizeSub*sizeof(double complex));

  matrix_trans(n, sizeSub, LREP->X1, transX, LREP->flag);
  matrix_trans(n, sizeSub, LREP->Y1, transY, LREP->flag);

  // caluclate denominator first
  // one norm of K&M already calculated upon setup
  for(int i = 0; i < sizeSub; i++) {
    denominator = CMPLX(LREP->oneNorm,0.0) + LREP->eValSort[i];

    // calculate nominator
    rsb_SPMV(LREP->K, &transX[i*n], work1, LREP->flag);

    copy_vec(n, &transY[i*n], work2, LREP->flag);
    zscale_vec(n, work2, LREP->eValSort[i]);
    sub_vec(n, work1, work2, LREP->flag);
    nominator = asum(n, work1, LREP->flag);

    rsb_SPMV(LREP->M, &transY[i*n], work1, LREP->flag);

    copy_vec(n, &transX[i*n], work2, LREP->flag);
    zscale_vec(n, work2, LREP->eValSort[i]);
    sub_vec(n, work1, work2, LREP->flag);
    nominator += asum(n, work1, LREP->flag);

    LREP -> resNorm[i] = cabs(nominator/denominator);
  }

  safe_free( work1 );
  safe_free( work2 );
  safe_free( transX );
  safe_free( transY );
  
}

void sp_get_residual_norm(struct sp_lrep_t *LREP) {
  switch(LREP->flag) {
  case 'Z':
    sp_get_residual_normz(LREP);
    break;

  default:
    sp_get_residual_normd(LREP);
    break;
  }
}
/*
  if(LREP->flag == 'D') {
    sp_get_residual_normd(LREP);
  } else if(LREP->flag == 'Z') {
    sp_get_residual_normz(LREP);
  }
  }*/


void sp_switch_eig_vec(struct sp_lrep_t *LREP) {
  // in the current iteration we saved the new eigenvectors in X1 and Y1 
  // however, in the next iteration those will be considered as
  // "the eigenvectors of the previous iteration"    
  // thus we move the current X1 and Y1 to X and Y and vice-versa   
  const int n = LREP->size;  
  const int sizeSub = LREP->sizeSub;   
   
  double *workd;
  double complex *workz;
  workd = xmalloc(n*sizeSub*sizeof(double));
  workz = xmalloc(n*sizeSub*sizeof(double complex));

  switch(LREP->flag) {
  case 'Z':
    copy_matrix(n, sizeSub, LREP->X, workz, LREP->flag);
    copy_matrix(n, sizeSub, LREP->X1, LREP->X, LREP->flag);
    copy_matrix(n, sizeSub, workz, LREP->X1, LREP->flag);

    copy_matrix(n, sizeSub, LREP->Y, workz, LREP->flag);
    copy_matrix(n, sizeSub, LREP->Y1, LREP->Y, LREP->flag);
    copy_matrix(n, sizeSub, workz, LREP->Y1, LREP->flag);

    safe_free( workz );
    safe_free( workd );
    break;

  default:
    copy_matrix(n, sizeSub, LREP->X, workd, LREP->flag);    
    copy_matrix(n, sizeSub, LREP->X1, LREP->X, LREP->flag);
    copy_matrix(n, sizeSub, workd, LREP->X1, LREP->flag);   
    
    copy_matrix(n, sizeSub, LREP->Y, workd, LREP->flag);    
    copy_matrix(n, sizeSub, LREP->Y1, LREP->Y, LREP->flag);
    copy_matrix(n, sizeSub, workd, LREP->Y1, LREP->flag);   
   
    safe_free( workd );
    safe_free( workz );
    break;
  }
}

/*
  if(LREP->flag == 'D') {
    copy_matrix(n, sizeSub, LREP->X, workd, LREP->flag);    
    copy_matrix(n, sizeSub, LREP->X1, LREP->X, LREP->flag);
    copy_matrix(n, sizeSub, workd, LREP->X1, LREP->flag);   
    
    copy_matrix(n, sizeSub, LREP->Y, workd, LREP->flag);    
    copy_matrix(n, sizeSub, LREP->Y1, LREP->Y, LREP->flag);
    copy_matrix(n, sizeSub, workd, LREP->Y1, LREP->flag);   
   
    safe_free( workd );
    safe_free( workz );

  } else if(LREP->flag == 'Z') {
    copy_matrix(n, sizeSub, LREP->X, workz, LREP->flag);
    copy_matrix(n, sizeSub, LREP->X1, LREP->X, LREP->flag);
    copy_matrix(n, sizeSub, workz, LREP->X1, LREP->flag);

    copy_matrix(n, sizeSub, LREP->Y, workz, LREP->flag);
    copy_matrix(n, sizeSub, LREP->Y1, LREP->Y, LREP->flag);
    copy_matrix(n, sizeSub, workz, LREP->Y1, LREP->flag);

    safe_free( workz );
    safe_free( workd );
  }
 
      
  }*/


void sp_lopb4dcg(sp_lrep_t *LREP) {

  int mult = 0;
  //  double s[3*LREP->sizeSub];
  int breaking = 0;


  printf("started LOPB4DCG.\n");
  sp_get_search_dir(LREP);

  printf("got search dir.\n");

  while(LREP->iter < LREP->maxIter) {
    mult = (0 == LREP->iter) ? 2 : 3;
    breaking = 0;
    
    sp_setup_P(LREP);
    sp_setup_Q(LREP);

    //    printf("pq\n");
    
    sp_apply_precond(LREP);
    //    printf("precond\n");

    sp_setup_U(LREP);
    sp_setup_V(LREP);
    //    printf("uv\n");

    get_q(LREP->size, mult*LREP->sizeSub, LREP->U, LREP->flag);
    get_q(LREP->size, mult*LREP->sizeSub, LREP->V, LREP->flag);
    //    printf("ortho\n");
    
    sp_setup_W(LREP);
    matrix_inverse(mult*LREP->sizeSub, LREP->W, LREP->flag);
    //calc_sing_val(mult*LREP->nev, mult*LREP->nev, s, LREP->W, LREP->flag);
    //sp_setup_W(LREP);
    //matrix_inverse(mult*LREP->nev, LREP->W, LREP->flag);
    sp_setup_Hsr(LREP);
    //    printf("Hsr\n");

    calc_eig_val(2*mult*LREP->sizeSub, LREP->Hsr, LREP->eValsr,
		  LREP->eVecsr, LREP->flag);
    //    printf("got eigenvalues.\n");

    sp_sort_eig(LREP);
    //    printf("sorted eigenvalues.\n");

    sp_split_eig_vec(LREP);

    sp_compute_eig_vec(LREP);
    sp_normalize_eig_vec(LREP);

    sp_get_residual_norm(LREP);

    printf("Eigenvalues at iteration %d:\n", LREP->iter+1);
    print_array_z(LREP->nev, 1, LREP->eValSort);
    printf("Residual Norm:\n");
    print_array_d(LREP->nev, 1, LREP->resNorm);
    printf("\n");

    sp_switch_eig_vec(LREP);

    // Assign eigenvalues as new search directions
#pragma omp parallel for
    for(int i = 0; i < LREP->sizeSub; i++)
      LREP->eigVal[i+i*LREP->sizeSub] = LREP->eValSort[i];

    LREP->iter += 1;

    breaking = 0;

    for(int i = 0; i < LREP->nev; i++)
      breaking += (LREP->resNorm[i] < LREP->tol) ? 1 : 0;
    /*
      if(LREP->resNorm[i] < LREP->tol)
	breaking += 1;
	}*/

    if(breaking == LREP->nev) 
      break;

  }

}
