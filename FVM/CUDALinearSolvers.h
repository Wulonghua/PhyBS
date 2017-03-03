#pragma once

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// CUDA Runtime
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_functions.h>

// Using updated (v2) interfaces for CUBLAS and CUSPARSE
#include <cusparse.h>
#include <cusolverSp.h>
#include <cublas_v2.h>

#include <helper_functions.h>  
#include <helper_cuda.h>

#include "GlobalHelpers.h"

class CUDALinearSolvers
{
public:
	CUDALinearSolvers(int rowNum, int nnz);
	~CUDALinearSolvers();

	void conjugateGradient(float *d_val, int *d_row, int *d_col, float *d_r, float *d_x);
	void directCholcuSolver(float *d_val, int *d_row, int *d_col, float *d_b, float *d_x);
	void chebyshevSemiIteration(float *d_val, int *d_row, int *d_col, int *d_diagonalIdx, float *d_b, float rho, float **d_y);
	void jacobiIteration(float *d_val, int *d_row, int *d_col, int *d_diagonalIdx, float *d_b, float **d_y);
	cublasHandle_t & getcuBlasHandle() { return cublasHandle; }




private:

	// split the matrix to diagonal and off-diagonal parts, A only remains the negative of off-diagonal parts,
	// B_1 stores the inverse of the diagonal parts.
	void splitCSRMatJacobi(float *Valptr, int m, int *diagonalIdx, float *B_1);

	int N;        // number of A's rows
	int nz;       // number of A's non-zero entries.

	/* CUBLAS context */
	cublasHandle_t cublasHandle;
	cublasStatus_t cublasStatus;

	/* Create CUSPARSE context */
	cusparseHandle_t cusparseHandle;
	cusparseStatus_t cusparseStatus;

	cusolverSpHandle_t cusolverHandle;
	cusolverStatus_t cusolverStatus;

	/* Description of the A matrix*/
	cusparseMatDescr_t m_descr;

	// followings are temped device memory for later conjugate gradient
	float *d_p;
	float *d_Ax;

	// followings are temped device memory for later Chebyshev Jacobi solver.
	float *d_y1;
	float *d_y2;

	// inverse of the diagonal of csrMat, used for jacobi iterative method.
	float *d_B_1;

	int m_N_threadsPerBlock;
	int m_N_blocksPerGrid;

};



