/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
/*
 *  Functions are modified from the NVIDIA official example code.
 *  Modified by Longhua Wu
 */
/*
 * This sample implements a preconditioned conjugate gradient solver on
 * the GPU using CUBLAS and CUSPARSE.  Relative to the conjugateGradient
 * SDK example, this demonstrates the use of cusparseScsrilu0() for
 * computing the incompute-LU preconditioner and cusparseScsrsv_solve()
 * for solving triangular systems.  Specifically, the preconditioned
 * conjugate gradient method with an incomplete LU preconditioner is
 * used to solve the Laplacian operator in 2D on a uniform mesh.
 *
 * Note that the code in this example and the specific matrices used here
 * were chosen to demonstrate the use of the CUSPARSE library as simply
 * and as clearly as possible.  This is not optimized code and the input
 * matrices have been chosen for simplicity rather than performance.
 * These should not be used either as a performance guide or for
 * benchmarking purposes.
 *
 */

#pragma once

#include "CUDALinearSolvers.h"

CUDALinearSolvers::CUDALinearSolvers(int rowNum, int nnz):
N(rowNum), nz(nnz)
{
	/* Create CUBLAS context */
    cublasHandle = 0;
    cublasStatus = cublasCreate(&cublasHandle);

    checkCudaErrors(cublasStatus);

    /* Create CUSPARSE context */
	cusparseHandle = 0;
    cusparseStatus = cusparseCreate(&cusparseHandle);

    checkCudaErrors(cusparseStatus);

    /* Description of the A matrix*/
    descr = 0;
    cusparseStatus = cusparseCreateMatDescr(&descr);

    checkCudaErrors(cusparseStatus);
	
	/* Define the properties of the matrix */
    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

	checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(float)));
	
}

CUDALinearSolvers::~CUDALinearSolvers()
{
	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);
	cudaFree(d_p);
	cudaFree(d_Ax);
}

void CUDALinearSolvers::conjugateGradient(float *d_val, int *d_row, int *d_col, float *d_r, float *d_x)
{
	float alpha = 1.0;
	float alpham1 = -1.0;
	float beta = 0.0;
	float r0 = 0.0;
	float a, b, na, r1, dot;
	const float tol = 1e-5f;
	const int max_iter = 10000;

	// compute A*x
	cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

	// compute d_r =  alpham1* d_Ax + d_r = d_r - d_Ax;  
	cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);

	// compute r1 = dot(d_1,d_1)
	cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

	int k = 1;

	while (r1 > tol*tol && k <= max_iter)
	{
		if (k > 1)
		{
			b = r1 / r0;
			cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
			cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
		}
		else
		{
			// d_p = d_r;
			cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
		}

		cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
		cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
		a = r1 / dot;

		cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
		na = -a;
		cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

		r0 = r1;
		cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		cudaThreadSynchronize();
		printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
		k++;
	}
	
}

