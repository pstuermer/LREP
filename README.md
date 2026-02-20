# LOPB4DCG-LREP
Implements a locally-optimal preconditioned block 4-dimensional conjugate gradient method to solve a linear response problem. Implementation based on https://epubs.siam.org/doi/10.1137/110838972 and https://link.springer.com/article/10.1007/s10543-014-0472-6.

This was my first scientific computing project and the predecessor to my LOBPCG implementation, which replaced this approach with a gradient-based method achieving 100x speedups.

## 1. What does this software do?
This software is a CPU-based Linear Response Eigenvalue Problem (LREP) Solver utilizing a conjugate-gradient type algorithm (https://epubs.siam.org/doi/10.1137/110838972 and https://link.springer.com/article/10.1007/s10543-014-0472-6). Usually, this kind of problem can be easily solved using diagonalization routines included in ARPACK or LAPACK. However, as this is the first piece of scientific software I wrote in C, I wanted to challenge myself and implement a non-standard algorithm. **This also means that the code itself is a learning project and contains a lot of beginner mistakes, and thereby is not recommended to be used as is. I plan on refactoring a majority of the code and extending it by the end of my PhD.**

For dense matrices, matrix-operations are implemented using OpenBlas (soon to be replaced with BLIS). For sparse matrices, they are implemented using the Recursive Sparse Blocks Format (https://librsb.sourceforge.net/). Currently represents operators via finite-difference on periodic boundary conditions using a plane-wave representation.

A Linear Response Eigenvalue Problem occurs if one were to allow the groundstate $\psi_0$ of a Density-Functional Theory type calculation to have quasiparticle excitations with energy $\hbar\omega$. The resulting Eigenvalue Problem can be written as 
```math
\begin{bmatrix}0&K\\M&0\end{bmatrix}\mathbf{u} =\hbar\omega\mathbf{u},
```
where both $K$ and $M$ are symmetric and at least one of them is positive semi-definite. In my field of ultra-cold Bose gases the 'Density-Functional Theory type calculation' is replaced by the so-called Gross-Pitaevskii equation or non-linear Schroedinger equation, assuming $T=0$.

Because of this, the code comes in two parts:
  - Construcing the matrices $K$ and $M$ in splrep.c/.h
  - Solving the resulting Eigenvalue Problem in splobp4dcg.c/.h

## 2. Works using LOPB4DCG-LREP
  - 'Mixed Bubbles in a one-dimensional Bose-Bose mixture': https://arxiv.org/abs/2207.11334

## 3. To-Do:
  - find a better name
  - write routine for dense matrices
  - replace OpenBlas by BLIS
  - replace standard dgemm routines by symmetric routines
  - parallelize kronecker product when constructing higher-dimensional differential operators
  - allow for remembering a solution of conjugate gradient preconditioner as initial guess for next iteration
  - allow for assigning function-pointers to assign different matrix-multiplication routines depending on some parameters
  - include more asserts
  - write unit-tests
  - replace if statements with function pointers at start of routine or switch statements
  - implement other finite difference schemes
  - implement LOBPCG for non-real $\psi_0$ (i.e. rotation in the case of a Bose gas)
  - write functions for reading and writing to file
  - allow for long-range interaction
  - structure factor
  - density of states
  - what about momentum space excitation spectrum (would immediately lead to non-symmetric $K$ and $M$)

