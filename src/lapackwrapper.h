#include "gen.h"
#include "openblaswrapper.h"
// *******************************************************************
// Includes wrapper routines for external lapack routines:                              
//                                                                                      
// - p_inverse:        Calculates the Pseudoinverse of a mxn 
//                     matrix A and save it in B
// - copy_matrix:      Copies a mxn matrix A to target                                    
// - get_1_norm:       Calculates the 1-norm of a mxn matrix                    
// - get_q:            Performs a QR-factorization of a mxn matrix
//                     A and then saves Q in A 
// - calc_sing_val:    Performs a singular value decomposition of 
//                     a mxn matrix A and saves         
//                     singular values in s. If needed vectors may 
//                     also be computed and saved         
// - calc_eig_val:     Calculates eigenvalues and righthand eigenvectors 
//                     of a nxn matrix    
//                     A and saves them in EigVal and EigVec respectively                   
// - matrix_inverse:   Calculates the inverse of a mxm matrix A via 
//                     LU-decomposition and      
//                     replaces the inverse with A
//      
// - Additional Note: The library is based on Fortran, which is column-major.           
//                    However, C is row-major. Thus we need to transpose                
//                    matrices before acting any routine on them.                       
//                                                                                      
// - Documentation about each library routine can be found on netlib.                   
//                                                                                      
// *******************************************************************
                                                                                        
/*---- IMPORT LAPACK FORTRAN ROUTINES ----*/                                            
extern void dgetrf_(int *M, int *N, double *A, int *lda,                                
                    int *IPIV, int *INFO);                                              
                                                                                        
extern void dgetri_(int *N, double *A, int *lad, int *IPIV,                             
                    double *WORK, int *lwork, int *INFO);                               
                                                                                        
extern void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A,                         
                   int *LDA, double *WR, double *WI, double *VL,                        
                   int *LDVL, double *VR, int *LDVR, double *WORK,                      
                   int *LWORK, int *INFO);                                              
                                                                                        
extern void dgesdd_(char *JOBZ, int *M, int *N, double *A, int *LDA,                    
                    double *S, double *U, int *LDU, double *VT,                         
                    int *LDVT, double *WORK, int *LWORK, int *IWORK,                    
                    int *INFO);                                                         
extern void dgeqrf_(int *M, int *N, double *A, int *LDA, double *TAU,                   
                    double *WORK, int *LWORK, int *INFO);                               
                                                                                        
extern void dorgqr_(int *M, int *N, int *K, double *A, int *LDA,                        
                    double *TAU, double *WORK, int *LWORK, int *INFO);                  
                                                                                        
extern void dormqr(char *SIDE, char *TRANS, int *M, int *N, int *K,                     
                   double *A, int *LDA, double *TAU, double *C,                         
                   int *LDC, double *WORK, int *LWORK, int *INFO);                      
                                                                                        
extern double dlange_(char *NORM, int *M, int *N, double *A,                            
                      int *LDA, double *WORK);                                          
                                                                                        
extern void dlacpy_(char *UPLO, int *M, int *N, double *A,                              
                    int *LDA, double *B, int *LDB);                                     
                                                                                        
extern void dgelss_(int *M, int *N, int *NRHS, double *A,                               
                    int *LDA, double *B, int *LDB, double *S,                           
                    double *RCOND, int *RANK, double *WORK,                             
                    int *LWORK, int *INFO);



extern void zgetrf_(int *M, int *N, double complex *A, int *lda,
		    int *IPIV, int *INFO);

extern void zgetri_(int *N, double complex *A, int *lad, int *IPIV,
		    double complex *WORK, int *lwork, int *INFO);

extern void zgeev_(char *JOBVL, char *JOBVR, int *N, double complex *A,
		   int *LDA, double complex *W, double complex *VL,
		   int *LDVL, double complex *VR, int *LDVR,
		   double complex *WORK, int *LWORK,
		   double *ROWRK, int *INFO);

extern void zgesdd_(char *JOBZ, int *M, int *N, double complex *A, int *LDA,
		    double *S, double complex *U, int *LDU,
		    double complex *VT, int *LDVT, double complex *WORK,
		    int *LWORK, double *RWORK, int *IWORK, int *INFO);

extern void zgeqrf_(int *M, int *N, double complex *A, int *LDA,
		    double complex *TAU, double complex *WORK,
		    int *LWORK, int *INFO);

extern void zungqr_(int *M, int *N, int *K, double complex *A, int *LDA,
		    double complex *TAU, double complex *WORK,
		    int *LWORK, int *INFO);
/*
extern void zormqr(char *SIDE, char *TRANS, int *M, int *N, int *K,
		   double *A, int *LDA, double *TAU, double *C,
		   int *LDC, double *WORK, int *LWORK, int *INFO);
*/
extern double zlange_(char *NORM, int *M, int *N, double complex *A,
		      int *LDA, double *WORK);

extern void zlacpy_(char *UPLO, int *M, int *N, double complex *A,
		    int *LDA, double complex *B, int *LDB);

extern void zgelss_(int *M, int *N, int *NRHS, double complex *A,
		    int *LDA, double complex *B, int *LDB, double complex *S,
		    double *RCOND, int *RANK, double complex *WORK,
		    int *LWORK, double complex *RWORK, int *INFO);


void *dp_inverse(const int m, const int n, double *matrix, double *target);
void *dcopy_matrix(const int m, const int n, double *matrix, double *target);    
double dget_1_norm(const int m, const int n, double *matrix);  
void *dget_q(const int m, const int n, double *matrix);     
void *dcalc_sing_val(const int m, const int n, double *sVal, double *matrix);   
void *dcalc_eig_val(const int n, double *matrix, double *eigVal, double *eigVec); 
void *dmatrix_inverse(const int m, double *matrix);


void *zp_inverse(const int m, const int n, double complex *matrix,
		 double complex *target);
void *zcopy_matrix(const int m, const int n, double complex *matrix,
		   double complex *target);
void *zget_q(const int m, const int n, double complex *matrix);
void *zcalc_sing_val(const int m, const int n, double *sVal, double complex *matrix);
void *zcalc_eig_val(const int n, double complex *matrix,
		    double complex *eigVal, double complex *eigVec);
void *zmatrix_inverse(const int m, double complex *matrix);


void *p_inverse(const int m, const int n, void *matrix,
		void *target, const char flag);
void *copy_matrix(const int m, const int n, void *matrix,
		  void *target, const char flag);
void *get_q(const int m, const int n, void *matrix,
	    const char flag);
void *calc_sing_val(const int m, const int n, double *sVal, void *matrix,
		    const char flag);
void *calc_eig_val(const int m, void *matrix, void *eigVal, void *eigVec,
		   const char flag);
void *matrix_inverse(const int m, void *matrix, const char flag);

double d_norm(const int m, const int n, double *matrix, char NORM);
double z_norm(const int m, const int n, double complex *matrix, char NORM);
double matrix_norm(const int m, const int n, void *matrix, char NORM,
		   const char flag);
