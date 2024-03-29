#pragma once

#include "gen.h"
#include "coo.h"
#include "grid.h"

// *******************************************************************
// Implements routines that are necessary to build the eigenvalue problem
// to be solved. This consists of some sparse differentialoperator in
// N-dimensions, plus some terms for trap geometries or interactions as
// detailed in the README.txt.
//
// Routines include:
// - setup_diff:    chooses one of the following three routines to build
//                  the differentiation operator based on the dim input
// - diff_1D:       builds a 1D spectral differentation operator
// - diff_2D:       builds a 2D spectral differentation operator
//                  using kronecker product
// - diff_3D:       builds a 3D spectral differentation operator 
//                  using kronecker product
// - sin_diff:      builds a dense spectral differentation operator
// - load_wf:       loads complex valued wavefunction based on input string
// - get_apsi:      calculates absolute valued wavefunction
// - cut_apsi:      cuts out a border of apsi and leaves the rest
// - setup_trap:    chooses one of the following three routines to build
//                  the diagonal trap coo_t matrix based on a grid, 
//                  dim and trapchoice input
// - trap_1D:       builds a 1D trap
// - trap_2D:       builds a 2D trap
// - trap_3D:       builds a 3D trap
// - setup_Uint:    based on dimensionality builds the diagonal 
//                  interaction term
// - d_Uint:        based on intchoice chooses the right interaction
// - setup_B:       based on dimensionality builds the second diagonal
//                  interaction term
// - b_int:         based on intchoice chooses the right interaction
// - setup_KM:      adds diff-op, trap and interactions together
//                  to build K & M matrices
// - init_eig_vec:  initializes eigenvectors within a range
// *******************************************************************

struct coo_t *setup_diff2(const double *ln, const int *N, const int dim,
			 const char flag);
struct coo_t *diff_1D(const double *ln, const int *N, const char flag);
struct coo_t *diff_2D(const double *ln, const int *N, const char flag);
struct coo_t *diff_3D(const double *ln, const int *N, const char flag);
void *fourier_diff2(const double lx, const int N, double *diffMatrix);
void fourier_diff1(const double lx, const int N, double *diffMatrix);


void *load_wf(double complex *wf, const int N, 
	      const int factor, char fileName[]);
void *load_wf_2C(double complex *wf1, double complex *wf2,
		 const int N, const int factor, char fileName[]);
double complex *get_psi_sq(double complex *wf, const int N);
double complex *get_ccpsi_sq(double complex *wf, const int N);
double *get_apsi(double complex *wf, const int N);
double *cut_apsi(double *apsi, const int *N, const int dim);
double complex *cut_wf(double complex *wf, const int *N, const int dim);


void setup_mu(double *diffMatrix, double *apsi1, double *apsi2, double mu[2],
	      const int intChoice, const double *intParam, grid_t *sys);
void *setup_trap(coo_t *cooTrap, const int *N, const double *param, grid_t *grid,
		 const int dim, const int trapChoice);
double trap_1D(const double *param, grid_t *sys, int *index, const int trapChoice);
double trap_2D(const double *param, grid_t *sys, int *index, const int trapChoice);
double trap_3D(const double *param, grid_t *sys, int *index, const int trapChoice);


void *setup_Uint(coo_t *cooUint, double *apsi1, double *apsi2, const int *N,
		 const int dim, const double mu, const double *intParam, 
		 const int intChoice, const int sysComp, int comp);
void *setup_Uint1C(coo_t *cooUint, double *apsi, const int *N, const int dim,
		   const double mu, const double *intParam, const int intChoice);
void* setup_Uint2C(coo_t *cooUint, double *apsi1, double *apsi2, const int *N,
		   const int dim, const double mu, const double *intParam,
		   const int intChoice, int comp);
double d_Uint1C(double apsi, const double mu, const double *intParam,
	     const int intChoice);
double d_Uint2C(double apsi1, double apsi2, const double mu,
		const double *intParam, const int intChoice, int comp);


void *setup_Bint(coo_t *coo_Bint, double *apsi1, double *apsi2,
		 const int *N, const int dim, const double *intParam,
		 const int intChoice, const int sysComp, int comp);
void *setup_Bint1C(coo_t *coo_Bint, double *apsi, const int *N,
		   const int dim, const double *intParam, const int intChoice);
void *setup_Bint2C(coo_t *coo_Bint, double *apsi1, double *apsi2,
		   const int *N, const int dim, const double *intParam, 
		   const int intChoice, int comp);
double d_Bint1C(double apsi, const double *intParam, const int intChoice);
double d_Bint2C(double apsi1, double apsi2, const double *intParam,
		const int intChoice, int comp);

void *setup_C(coo_t *cooC, double *apsi1, double *apsi2, const int *N,
	      const int dim, const double *intParam, const int intChoice);
double d_Cint(double dapsi, double apsi1, double apsi2,
	      const double *intParam, const int intChoice);


void *setup_KM1C(coo_t *K, coo_t *diff, coo_t *trap, coo_t *uint,
	       coo_t *B1, const char mode, const int nnz,
	       const char flag);
void* setup_KM2C(coo_t *KM, coo_t *diff, coo_t *trap, coo_t *Uint1, coo_t *Uint2,
		 coo_t *B1, coo_t *B2, coo_t *C, const char mode, const int nnz,
		 const char flag);

void *init_eig_vecd(const int N, const int nev, double *vec,
		    double low, double high);
void *init_eig_vecz(const int N, const int nev, double complex *vec,
		    double low, double high);
void *init_eig_vec(const int N, const int nev, void *vec,
		   double low, double high, const char flag);

void add_mode2C(const int dim, const int *N, double complex *mode,
	      const char wfString[], const char saveString[]);
