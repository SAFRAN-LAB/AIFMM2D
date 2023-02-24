This library solves linear systems Ax=b, arising out of N-body problems in 2D. Such matrices A arising out of N-body problems are rank-structured. By exploiting this rank-structure, a fast direct sovler is built, that is termed Algebraic Inverse Fast Multipole Method (AIFMM).

For more details please refer[[1](https://arxiv.org/pdf/2301.12704.pdf)].

The location of particles is to be defined in the userkernel class in the vector "gridPoints". Currently particles are distributed uniformly in the square [-L,L]^{2}.

The matrix entries are to be defined in the function "getMatrixEntry(i,j)" of the kernel.hpp file.
 Currently, 1/r and Hankel function are defined. Use DTYPE_FLAG in Makefile2D.mk to choose between 1/r and Hankel.

The rhs vector is to be defined in VectorXd "true_Ax" of the testAIFMM.cpp file. Currently for test purposes, it is defined as a vector containing ones and zeros at random locations.
