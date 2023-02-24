This library solves linear systems Ax=b, arising out of N-body problems in 2D. Such matrices A arising out of N-body problems are rank-structured. By exploiting this rank-structure, a fast direct sovler is built, that is termed Algebraic Inverse Fast Multipole Method (AIFMM).

For more details please refer[[1](https://arxiv.org/pdf/2301.12704.pdf)].

The location of particles is to be defined in the userkernel class in the vector "gridPoints". Currently particles are distributed uniformly in the square [-L,L]^{2}.

The matrix entries are to be defined in the function "getMatrixEntry(i,j)" of the kernel.hpp file.
 Currently, 1/r and Hankel function are defined. Use DTYPE_FLAG in Makefile2D.mk to choose between 1/r and Hankel.

The rhs vector is to be defined in VectorXd "true_Ax" of the testAIFMM.cpp file. Currently for test purposes, it is defined as a vector containing ones and zeros at random locations.

The testAIFMM.cpp takes the following inputs at run time:

sqrtRootN: square root of N.

L: half side length of the computational square domain centered at 0

TOL_POW: the compression tolerance used by ACA routine

nParticlesInLeafAlong1D: maximum number of particles in leaf in a single dimension, i.e., in the process of creating hierarchical tree, a box is considered to be a leaf if it contains atmost nParticlesInLeafAlong1D*nParticlesInLeafAlong1D particles.

A sample output looks like this:

N: 10000

nLevelsUniform: 5

Max rank of compressible sub-matrices: 29

Assemble time: 1.79592

Factorization time: 11.9208

Solve time: 0.11616

Relative forward error in the solution: 1.68157e-05
