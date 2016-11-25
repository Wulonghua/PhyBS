#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// CUDA Runtime
#include <cuda_runtime.h>

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
		cublasHandle_t & getcuBlasHandle() { return cublasHandle; }




private:
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
    cusparseMatDescr_t descr;

	// followings are temped device memory for later conjugate gradient solver.
	float *d_p;
	float *d_Ax;

};



