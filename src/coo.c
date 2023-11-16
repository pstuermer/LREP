
#include "coo.h"

struct coo_t *coo_malloc(const int nonz, const int numRows, const int numCols,
			   const char flag) {
  coo_t *matrix = xmalloc(sizeof(struct coo_t));

  matrix -> nnz = nonz;
  matrix -> rows = numRows;
  matrix -> cols = numCols;
  matrix -> last = 0;
  matrix -> flag = flag;

  if(matrix->flag == 'D') {
    matrix->dval = xmalloc(nonz*sizeof(double));
    matrix->zval = xmalloc(0*sizeof(double complex));
  } else if (matrix->flag == 'Z') {
    //    printf("%d\n",nonz);
    matrix->dval = xmalloc(0*sizeof(double));
    matrix->zval = xmalloc(nonz*sizeof(double complex));
  }
  //  matrix -> val = xmalloc(nonz*sizeof(double));
  matrix -> col = xmalloc(nonz*sizeof(int));
  matrix -> row = xmalloc(nonz*sizeof(int));
  matrix -> numElRow = calloc(numRows, sizeof(int));

  return matrix;
}

void coo_free(coo_t *matrix) {
  assert(matrix != NULL);
  
  safe_free( matrix -> dval );
  safe_free( matrix -> zval );
  //  safe_free( matrix -> val );
  safe_free( matrix -> col );
  safe_free( matrix -> row );
  safe_free( matrix -> numElRow );

  safe_free( matrix );

  matrix = NULL;
}

void coo_reset(coo_t *matrix) {
  matrix -> last = 0;
}

void coo_insert(coo_t *matrix, double real, double imag, int colPos, int rowPos) {
  // requires that the input to be inserted comes after in the matrix
  // than all previous inputs (i.e. sorted input)

  assert(matrix -> last < matrix->nnz);

  double complex val = CMPLX(real, imag);
  if(matrix->flag == 'D') {
    matrix->dval[matrix->last] = real;
  } else if (matrix->flag == 'Z') {
    matrix->zval[matrix->last] = val;
  }

  matrix -> col[matrix->last] = colPos;
  matrix -> row[matrix->last] = rowPos;
  matrix -> last += 1;
}

void coo_get_num_row(coo_t *matrix) {
  int currRow = 0;
  for(int i = 0; i < matrix->nnz; i++) {
    if (matrix->row[i] == currRow) {
      matrix -> numElRow[currRow] += 1;
    } else {
      currRow = matrix -> row[i];
      matrix -> numElRow[currRow] += 1;
    }
  }
}

void dcoo_kron(coo_t *matrixLeft, coo_t *matrixRight, coo_t *result) { 
  int sumElLeft, sumElRight;
  int posKron;
  int rowPos, colPos;
  int colOffset;
  double leftVal, val;
  posKron = sumElLeft = 0;

  if(matrixLeft->numElRow[0] == 0 || 
     matrixLeft->numElRow[0] > matrixLeft->nnz)
    coo_get_num_row(matrixLeft);
  if(matrixRight->numElRow[0] == 0 ||
     matrixRight->numElRow[0] > matrixRight->nnz)
    coo_get_num_row(matrixRight);
                                                                   
  for(int i = 0; i < matrixLeft->rows; i++) { 
    sumElRight = 0;
    for(int k = 0; k < matrixRight->rows; k++) {
      rowPos = i*matrixLeft->rows+k;    
      for(int j = 0; j < matrixLeft->numElRow[i]; j++) {     
        leftVal = matrixLeft->dval[sumElLeft+j];     
        colOffset = matrixLeft->col[sumElLeft+j]*matrixLeft->cols;    
        for(int l = 0; l < matrixRight->numElRow[k]; l++) {    
          colPos = colOffset + matrixRight->col[sumElRight+l];   
          val = leftVal * matrixRight->dval[sumElRight+l];                                           
          coo_insert(result, val, 0, colPos, rowPos); 
          posKron += 1;      
        }      
      }                    
      sumElRight += matrixRight->numElRow[k];     
    }      
    sumElLeft += matrixLeft->numElRow[i];      
  }                  
} 


void dcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR) {  
  // The idea for addition of two COO matrices is simple
  // I have a result matrix as input, which has at most nnza+nnzb nonzero elements   
  // starting for both left and right at the first position, I incrementally 
  // increase the index indivudally for both on the following  
  // 1. if both have an element at the same position, we add the two values, 
  // check that the result is not zero, assign the position and increment
  // both a and b  
  // if a->row < b->row || a->row == b->row && a->col < b->col, we insert the value 
  // from a at the position of a, and increment a. The same for b.
  // If we reach the end of a or b, we exit the loop and assign the   
  // remaining values of the other matrix to the result matrix 
                                    
  int aPos, bPos;  
  aPos = 0;
  bPos = 0;
  double value;
  while( aPos < matrixA->nnz && bPos < matrixB->nnz) {
    // if both aPos and bPos give an element at the same position
    if((matrixA->col[aPos] == matrixB->col[bPos]) &&
       (matrixA->row[aPos] == matrixB->row[bPos])) {
      value = matrixA->dval[aPos] + matrixB->dval[bPos];

      if (value != 0) {
        coo_insert(matrixR, value, 0, matrixA->col[aPos], matrixA->row[aPos]);
        bPos++;         
        aPos++; 
      } else {    
        bPos++;  
        aPos++;  
      }      
      // if element in A is "earlier" than element in b 
    } else if ( (matrixA->row[aPos] < matrixB->row[bPos]) ||
                (matrixA->row[aPos] == matrixB->row[bPos] && 
                 matrixA->col[aPos] < matrixB->col[bPos]) ) {    
      coo_insert(matrixR, matrixA->dval[aPos], 0, 
		 matrixA->col[aPos], matrixA->row[aPos]);
      aPos++;     
    } else if ( (matrixA->row[aPos] > matrixB->row[bPos]) ||  
                (matrixA->row[aPos] == matrixB->row[bPos] &&   
                 matrixA->col[aPos] > matrixB->col[bPos]) ) { 
      coo_insert(matrixR, matrixB->dval[bPos], 0,
		 matrixB->col[bPos], matrixB->row[bPos]);   
      bPos++;                        
    }                    
  }
  // now one of the conditions has triggered, so we need to
  // add the remaining values of either a or b
  while (aPos < matrixA->nnz) {   
    coo_insert(matrixR, matrixA->dval[aPos], 0,
	       matrixA->col[aPos], matrixA->row[aPos]);
    aPos++;         
  }          
  while (bPos < matrixB->nnz) {    
    coo_insert(matrixR, matrixB->dval[bPos], 0,
	       matrixB->col[bPos], matrixB->row[bPos]); 
    bPos++;       
  }                                 
  matrixR->nnz = matrixR->last;      
}


void dcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR) {   
  int aPos, bPos;   
  aPos = 0;     
  bPos = 0;   
  double value;   
  while( aPos < matrixA->nnz && bPos < matrixB->nnz) {  
    // if both aPos and bPos give an element at the same position  
    if((matrixA->col[aPos] == matrixB->col[bPos]) &&
       (matrixA->row[aPos] == matrixB->row[bPos])) {  
      value = matrixA->dval[aPos] - matrixB->dval[bPos];
      if (value != 0) {    
        coo_insert(matrixR, value, 0, matrixA->col[aPos], matrixA->row[aPos]);   
        bPos++;           
        aPos++;      
      } else {        
        bPos++;    
        aPos++;       
      }      
      // if element in A is "earlier" than element in b  
    } else if ( (matrixA->row[aPos] < matrixB->row[bPos]) ||    
                (matrixA->row[aPos] == matrixB->row[bPos] &&      
                 matrixA->col[aPos] < matrixB->col[bPos]) ) {    
      coo_insert(matrixR, matrixA->dval[aPos], 0,
		 matrixA->col[aPos], matrixA->row[aPos]);
      aPos++;   
    } else if ( (matrixA->row[aPos] > matrixB->row[bPos]) || 
                (matrixA->row[aPos] == matrixB->row[bPos] && 
                 matrixA->col[aPos] > matrixB->col[bPos]) ) { 
      coo_insert(matrixR, -matrixB->dval[bPos], 0,
		 matrixB->col[bPos], matrixB->row[bPos]);   
      bPos++;              
    }       
  }               
  // now one of the conditions has triggered, so we need to     
  // add the remaining values of either a or b
  while (aPos < matrixA->nnz) {      
    coo_insert(matrixR, matrixA->dval[aPos], 0,
	       matrixA->col[aPos], matrixA->row[aPos]); 
    aPos++; 
  }                   
  while (bPos < matrixB->nnz) {      
    coo_insert(matrixR, -matrixB->dval[bPos], 0,
	       matrixB->col[bPos], matrixB->row[bPos]);   
    bPos++;          
  }                 
  matrixR->nnz = matrixR->last; 
}    
/*                                                        
double dcoo_1norm(coo_t *matrix) {  
  double *work;   
  double oneNorm;  
  oneNorm = 0; 
  work = xmalloc(matrix->cols*sizeof(double));
                                              
  for(int i = 0; i < matrix->nnz; i++) {      
    work[matrix->col[i]] += abs(matrix->dval[i]);
  }  
     
  for(int i = 0; i < matrix->cols; i++) {  
    oneNorm = MAX(oneNorm, work[i]);     
  }
                                                                                             
  return oneNorm; 
}
*/
/*
struct zcoo_t *zcoo_malloc(const int nonz, const int numRows, const int numCols) {
  zcoo_t *matrix = xmalloc(sizeof(zcoo_t));

  matrix -> nnz = nonz;
  matrix -> rows = numRows;
  matrix -> cols = numCols;
  matrix -> last = 0;

  matrix -> val = xmalloc(nonz*sizeof(double complex));
  matrix -> col = xmalloc(nonz*sizeof(int));
  matrix -> row = xmalloc(nonz*sizeof(int));
  matrix -> numElRow = xcalloc(numRows, sizeof(int));

  return matrix;
}

void zcoo_free(zcoo_t *matrix) {
  assert(matrix != NULL);

  safe_free( matrix -> val );
  safe_free( matrix -> col );
  safe_free( matrix -> row );
  safe_free( matrix -> numElRow );

  safe_free( matrix );
}

void zcoo_insert(zcoo_t *matrix, double complex val,
		 int colPos, int rowPos) {
  assert(matrix -> last < matrix-> nnz);

  matrix -> val[matrix->last] = val;
  matrix -> col[matrix->last] = colPos;
  matrix -> row[matrix->last] = rowPos;
  matrix -> last += 1;
}

void zcoo_get_num_row(zcoo_t *matrix) {
  int currRow = 0;
  for(int i = 0; i < matrix->nnz; i++) {
    if (matrix->row[i] == currRow) {
      matrix -> numElRow[currRow] +=1;
    } else {
      currRow = matrix -> row[i];
      matrix -> numElRow[currRow] += 1;
    }
  }
}
*/
void zcoo_kron(coo_t *matrixLeft, coo_t *matrixRight, coo_t *result) {
  int sumElLeft, sumElRight;
  int posKron;
  int rowPos, colPos;
  int colOffset;
  double complex leftVal, val;
  posKron = sumElLeft = 0;

  if(matrixLeft->numElRow[0] == 0 ||
     matrixLeft->numElRow[0] > matrixLeft->nnz)
    coo_get_num_row(matrixLeft);
  if(matrixRight->numElRow[0] == 0 ||
     matrixRight->numElRow[0] > matrixRight->nnz)
    coo_get_num_row(matrixRight);

  for(int i = 0; i < matrixLeft->rows; i++) {
    sumElRight = 0;
    for(int k = 0; k < matrixRight->rows; k++) {
      rowPos = i*matrixLeft->rows+k;
      for(int j = 0; j < matrixLeft->numElRow[i]; j++) {
	leftVal = matrixLeft->zval[sumElLeft+j];
	colOffset = matrixLeft->col[sumElLeft+j]*matrixLeft->cols;
	for(int l = 0; l < matrixRight->numElRow[k]; l++) {
	  colPos = colOffset + matrixRight->col[sumElRight+l];
	  val = leftVal * matrixRight->zval[sumElRight+l];
	  coo_insert(result, creal(val), cimag(val), colPos, rowPos);
	  posKron += 1;
	}
      }
      sumElRight += matrixRight->numElRow[k];
    }
    sumElLeft += matrixLeft->numElRow[i];
  }
}

void zcoo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR){
  // The idea for addition of two COO matrices is simple
  // I have a result matrix as input, which has at most nnza+nnzb nonzero elements
  // starting for both left and right at the first position, I incrementally
  // increase the index indivudally for both on the following
  // 1. if both have an element at the same position, we add the two values,
  // check that the result is not zero, assign the position and increment
  // both a and b
  // if a->row < b->row || a->row == b->row && a->col < b->col, we insert the value
  // from a at the position of a, and increment a. The same for b.
  // If we reach the end of a or b, we exit the loop and assign the
  // remaining values of the other matrix to the result matrix

  int aPos, bPos;
  aPos = 0;
  bPos = 0;
  double complex value;
  //  printf("ABR %d\t%d\t%d\n",matrixA->nnz,matrixB->nnz,matrixR->nnz);
  while( aPos < matrixA->nnz && bPos < matrixB->nnz) {
    // if both aPos and bPos give an element at the same position
    if((matrixA->col[aPos] == matrixB->col[bPos]) &&
       (matrixA->row[aPos] == matrixB->row[bPos])) {
      value = matrixA->zval[aPos] + matrixB->zval[bPos];
      if (value != 0) {
	coo_insert(matrixR, creal(value), cimag(value),
		   matrixA->col[aPos], matrixA->row[aPos]);
	bPos++;
	aPos++;
      } else {
	bPos++;
	aPos++;
      }
      //      printf("both equal\t");
      // if element in A is "earlier" than element in b
    } else if ( (matrixA->row[aPos] < matrixB->row[bPos]) ||
		(matrixA->row[aPos] == matrixB->row[bPos] &&
		 matrixA->col[aPos] < matrixB->col[bPos]) ) {
      coo_insert(matrixR, creal(matrixA->zval[aPos]),
		 cimag(matrixA->zval[aPos]),
		 matrixA->col[aPos], matrixA->row[aPos]);
      aPos++;
    } else if ( (matrixA->row[aPos] > matrixB->row[bPos]) ||
		(matrixA->row[aPos] == matrixB->row[bPos] &&
		 matrixA->col[aPos] > matrixB->col[bPos]) ) {
      coo_insert(matrixR, creal(matrixB->zval[bPos]),
		 cimag(matrixB->zval[bPos]),
		 matrixB->col[bPos], matrixB->row[bPos]);
      bPos++;
    }
    //    printf("%d\t",aPos+bPos);
  }
  // now one of the conditions has triggered, so we need to
  // add the remaining values of either a or b
  while (aPos < matrixA->nnz) {
    //    printf("a\t");
    coo_insert(matrixR, creal(matrixA->zval[aPos]),
	       cimag(matrixA->zval[aPos]),
	       matrixA->col[aPos], matrixA->row[aPos]);
    aPos++;
    //    printf("%d\t",aPos+bPos);
  }
  //  printf("\n");
  //  printf("apos %d\n",aPos);
  while (bPos < matrixB->nnz) {
    //    printf("b\t");
    coo_insert(matrixR, creal(matrixB->zval[bPos]),
	       cimag(matrixA->zval[aPos]),
	       matrixB->col[bPos], matrixB->row[bPos]);
    bPos++;
    //    printf("%d\t",aPos+bPos);
  }
  //  printf("\n");
  //  printf("bpos %d\n",bPos);
  //  printf("\n");
  matrixR->nnz = matrixR->last;
}

void zcoo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR) {
  int aPos, bPos;
  aPos = 0;
  bPos = 0;
  double complex value;
  while( aPos < matrixA->nnz && bPos < matrixB->nnz) {
    // if both aPos and bPos give an element at the same position
    if((matrixA->col[aPos] == matrixB->col[bPos]) &&
       (matrixA->row[aPos] == matrixB->row[bPos])) {
      value = matrixA->zval[aPos] - matrixB->zval[bPos];
      if (value != 0) {
	coo_insert(matrixR, creal(value), cimag(value),
		   matrixA->col[aPos], matrixA->row[aPos]);
	bPos++;
	aPos++;
      } else {
	bPos++;
	aPos++;
      }
      // if element in A is "earlier" than element in b
    } else if ( (matrixA->row[aPos] < matrixB->row[bPos]) ||
		(matrixA->row[aPos] == matrixB->row[bPos] &&
		 matrixA->col[aPos] < matrixB->col[bPos]) ) {
      coo_insert(matrixR, creal(matrixA->zval[aPos]),
		 cimag(matrixA->zval[aPos]),
		 matrixA->col[aPos], matrixA->row[aPos]);
      aPos++;
    } else if ( (matrixA->row[aPos] > matrixB->row[bPos]) ||
		(matrixA->row[aPos] == matrixB->row[bPos] &&
		 matrixA->col[aPos] > matrixB->col[bPos]) ) {
      coo_insert(matrixR, creal(-matrixB->zval[bPos]),
		 cimag(-matrixB->zval[bPos]),
		 matrixB->col[bPos], matrixB->row[bPos]);
      bPos++;
    }
  }
  // now one of the conditions has triggered, so we need to
  // add the remaining values of either a or b
  while (aPos < matrixA->nnz) {
    coo_insert(matrixR, creal(matrixA->zval[aPos]),
	       cimag(matrixA->zval[aPos]),
	       matrixA->col[aPos], matrixA->row[aPos]);
    aPos++;
  }
  while (bPos < matrixB->nnz) {
    coo_insert(matrixR, creal(-matrixB->zval[bPos]),
	       cimag(matrixB->zval[bPos]),
	       matrixB->col[bPos], matrixB->row[bPos]);
    bPos++;
  }
  matrixR->nnz = matrixR->last;
}

double coo_1norm(coo_t *matrix) {
  double *work;
  double oneNorm;
  oneNorm = 0;
  work = xmalloc(matrix->cols*sizeof(double));

  for(int i = 0; i < matrix->nnz; i++) {
    if(matrix->flag == 'D') {
      work[matrix->col[i]] += fabs(matrix->dval[i]);
    } else if(matrix->flag == 'Z') {
      work[matrix->col[i]] += cabs(matrix->zval[i]);
    }
  }

  for(int i = 0; i < matrix->cols; i++) {
    oneNorm = MAX(oneNorm, work[i]);
  }

  return oneNorm;
}

void coo_kron(coo_t *left, coo_t *right, coo_t *result) {
  if(left->flag != right->flag ||
     left->flag != result->flag ||
     right->flag != result->flag) {
    printf("Datatypes in coo_kron don't fit.\n");
    exit( 1 );
  }

  if(left->flag == 'D') {
    dcoo_kron(left, right, result);
  } else if(left->flag == 'Z') {
    zcoo_kron(left, right, result);
  }
}

void coo_add(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR) {
  if(matrixA->flag != matrixB->flag ||
     matrixA->flag != matrixR->flag ||
     matrixB->flag != matrixR->flag) {
    printf("Datatypes in coo_add don't fit.\n");
    exit( 1 );
  }

  if(matrixA->flag == 'D') {
    dcoo_add(matrixA, matrixB, matrixR);
  } else if(matrixA->flag == 'Z') {
    zcoo_add(matrixA, matrixB, matrixR);
  }
}

void coo_sub(coo_t *matrixA, coo_t *matrixB, coo_t *matrixR) {
  if(matrixA->flag != matrixB->flag ||
     matrixA->flag != matrixR->flag ||
     matrixB->flag != matrixR->flag) {
    printf("Datatypes in coo_sub don't fit.\n");
    exit( 1 );
  }

  if(matrixA->flag == 'D') {
    dcoo_sub(matrixA, matrixB, matrixR);
  } else if(matrixA->flag == 'Z') {
    zcoo_sub(matrixA, matrixB, matrixR);
  }
}
