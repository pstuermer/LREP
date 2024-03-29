#pragma once

#include "gen.h"
#include "splrep.h"
#include "rsbwrapper.h"
#include "openblaswrapper.h"
#include "lapackwrapper.h"
#include "precond.h"

// *******************************************************************
// Implements all functions to solve a sparse linear response eigenvalue
// problem via locally optimal preconditioned block 4-dimensional 
// conjugate gradient:
// sp_get_search_dir:    only on first iteration, calculate eigenvalues
//                       of randomly created eigenvectors
// sp_setup_P:           calculate gradient search-direction for new X
// sp_setup_Q:           calculate gradient search-direction for new Y
// sp_apply_precond:     Apply preconditioner to new P and Q
// sp_setup_U:           setup U from X, X1 and P
// sp_setup_V:           setup V from Y, Y1 and Q
// sp_wetup_W:           calculates W from U and V, which is used for
//                       preserving symmetry in Hsr
// sp_setup_Hsr:         sets up the reduced structure preserving 
//                       hamiltonian of the linear response
//                       eigenvalue problem
// sp_sort_eig:          sort eigenvalues & -vectors and then only choose
//                       the k-smallest positive ones
// sp_split_eig_vec:     eigenvectors calculated are X and Y in one,
//                       so need to split them into their components
// sp_compute_eig_vec:   calculates new eigenvectors from reduced eigvec
// sp_normalize_eig_vec: Normalizes calculated eigenvectors to unity
// sp_get_residual_norm: calculates residual norm of new eigenvectors and
//                       values
// sp_switch_eig_vec:    switch current X&Y to X1&Y1, switch new eigen-
//                       vectors to X&Y
// *******************************************************************

void *sp_get_search_dird(struct sp_lrep_t *LREP);
void *sp_get_search_dirz(struct sp_lrep_t *LREP);
void *sp_get_search_dir(struct sp_lrep_t *LREP);

void *sp_setup_P(struct sp_lrep_t *LREP);
void *sp_setup_Q(struct sp_lrep_t *LREP);

void *sp_apply_precondd(struct sp_lrep_t *LREP);
void *sp_apply_precondz(struct sp_lrep_t *LREP);
void *sp_apply_precond(struct sp_lrep_t *LREP);

void *sp_setup_U(struct sp_lrep_t *LREP);
void *sp_setup_V(struct sp_lrep_t *LREP);
void *sp_setup_W(struct sp_lrep_t *LREP);

void *sp_setup_Hsrd(struct sp_lrep_t *LREP);
void *sp_setup_Hsrz(struct sp_lrep_t *LREP);
void *sp_setup_Hsr(struct sp_lrep_t *LREP);

void *sp_sort_eigd(struct sp_lrep_t *LREP);
void *sp_sort_eigz(struct sp_lrep_t *LREP);
void *sp_sort_eig(struct sp_lrep_t *LREP);

void *sp_split_eig_vec(struct sp_lrep_t *LREP);
void *sp_compute_eig_vec(struct sp_lrep_t *LREP);
void *sp_normalize_eig_vec(struct sp_lrep_t *LREP);

void *sp_get_residual_normd(struct sp_lrep_t *LREP);
void *sp_get_residual_normz(struct sp_lrep_t *LREP);
void *sp_get_residual_norm(struct sp_lrep_t *LREP);

void *sp_switch_eig_vec(struct sp_lrep_t *LREP);

void *sp_lopb4dcg(struct sp_lrep_t *LREP);


