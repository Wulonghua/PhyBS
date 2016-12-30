#pragma once

#include "CUDALinearSolvers.h"

// c_i = a_i * b_i
__global__ void multiplyArrayElementwise(float *a, float *b, int n, float *c)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	c[i] = a[i] * b[i];
}

// y2 = w2 * ( B_1 *(p + b) - y) + y
// p = C * y1
// see equation(12) in [Wang 15].
__global__ void updateChebyshevSemiIteration(float w, float *B_1, float *p, float *b, float *y, int n, float *y2)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	y2[i] = w *(B_1[i] * (p[i] + b[i]) - y[i]) + y[i];
}

CUDALinearSolvers::CUDALinearSolvers(int rowNum, int nnz):
N(rowNum), nz(nnz)
{
	/* Create CUBLAS context */
    //cublasHandle = 0;
    cublasStatus = cublasCreate(&cublasHandle);

    checkCudaErrors(cublasStatus);

    /* Create CUSPARSE context */
	//cusparseHandle = 0;
    cusparseStatus = cusparseCreate(&cusparseHandle);

    checkCudaErrors(cusparseStatus);

	/* Create CUSOLVER context*/
	cusolverStatus = cusolverSpCreate(&cusolverHandle);
	checkCudaErrors(cusolverStatus);

    /* Description of the A matrix*/
    //m_descr = 0;
	cusparseStatus = cusparseCreateMatDescr(&m_descr);

    checkCudaErrors(cusparseStatus);
	
	/* Define the properties of the matrix */
	cusparseSetMatType(m_descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(m_descr, CUSPARSE_INDEX_BASE_ZERO);

	checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_y1, N*sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_y2, N*sizeof(float)));

	m_N_threadsPerBlock = 64;
	m_N_blocksPerGrid = (N + m_N_threadsPerBlock - 1) / m_N_threadsPerBlock;

}

CUDALinearSolvers::~CUDALinearSolvers()
{
	cusparseDestroy(cusparseHandle);
	cublasDestroy(cublasHandle);
	cusolverSpDestroy(cusolverHandle);
	cudaFree(d_p);
	cudaFree(d_Ax);
	cudaFree(d_y1);
	cudaFree(d_y2);
}

void CUDALinearSolvers::directCholcuSolver(float *d_val, int *d_row, int *d_col, float *d_b, float *d_x)
{
	int reorder = 0; //no reordering
	int singularity = 0;
	cusolverStatus = cusolverSpScsrlsvchol(cusolverHandle, N, nz, m_descr, d_val, d_row, d_col, d_b, 1e-3, reorder, d_x, &singularity);
	checkCudaErrors(cusolverStatus);
}


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
void CUDALinearSolvers::conjugateGradient(float *d_val, int *d_row, int *d_col, float *d_r, float *d_x)
{
	float alpha = 1.0;
	float alpham1 = -1.0;
	float beta = 0.0;
	float r0 = 0.0;
	float a, b, na, r1, dot;
	const float tol = 1e-6f;
	const int max_iter = 10000;

	// compute A*x
	cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, m_descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

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

		cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, m_descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
		cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
		a = r1 / dot;

		cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
		na = -a;
		cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

		r0 = r1;
		cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
		cudaThreadSynchronize();
		//printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
		k++;
	}
	
}

void CUDALinearSolvers::chebyshevSemiIterativeSolver(float *d_val, int *d_row, int *d_col, float *d_B_1, float *d_b, float rho, float **d_y)
{
	//initialize y, y1
	// y = 0;
	// y1 = B_1 * (C*y+b)
	cudaMemset(*d_y, 0, N*sizeof(float));
	multiplyArrayElementwise<<<m_N_blocksPerGrid,m_N_threadsPerBlock>>>(d_B_1, d_b, N, d_y1);

	float omega = 2.0 / (2.0 - rho * rho);

	float alpha = 1.0;
	float beta = 0.0;
	//iteration 
	for (int k = 0; k < 450;++k)
	{
		cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, m_descr, d_val, d_row, d_col, d_y1, &beta, d_p);
		updateChebyshevSemiIteration<<<m_N_blocksPerGrid, m_N_threadsPerBlock>>>(omega, d_B_1, d_p, d_b, d_y, N, d_y2);
		cudaDeviceSynchronize();
		swap(d_y, &d_y1);
		swap(&d_y1, &d_y2);
		omega = 4.0 / (4.0 - rho*rho*omega);
	}

	swap(d_y,&d_y1);
}