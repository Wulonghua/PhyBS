#pragma once

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

class CUDALinearSolvers
{
public:
	CUDALinearSolvers(int rowNum, int nnz);
	~CUDALinearSolvers();

	void conjugateGradient(float *d_val, int *d_row, int *d_col, float *d_r, float *d_x);
	void directCholcuSolver(float *d_val, int *d_row, int *d_col, float *d_b, float *d_x);
	void chebyshevSemiIterativeSolver(float *d_val, int *d_row, int *d_col, float *d_B_1, float *d_b, float rho, float **d_y);
	cublasHandle_t & getcuBlasHandle() { return cublasHandle; }




private:
	void swap(float **a, float **b)
	{
		float *t = *a;
		*a = *b;
		*b = t;
	}

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

	int m_N_threadsPerBlock;
	int m_N_blocksPerGrid;

};



