#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include "CUDA_MATH.h"       // helper matrix operations 
#include "helper_math.h" // helper vector(float3) operations


class CUDAInterface
{
public:
	CUDAInterface(int num_nodes, const float *nodes, int num_tets, const int *tets, 
		const float youngs, const float nu, const float density, 
		int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval);
	~CUDAInterface();


private:
	  int n_nodes;
	  int n_tets;
	  int *d_tets;
	  float	*d_nodes;
	  float *d_masses;
	  float *d_ANs;
	  float *d_Dm_inverses;
	  float *d_mus;
	  float *d_lambdas;
	  float *d_Fhats;
	  float *d_Us;
	  float *d_Vs;

	  struct CSRmatrix{
		  int nnz;   // number of non-zero elements
		  int m;     // rows' number
		  float *d_Valptr;
		  int *d_Rowptr;
		  int *d_Colptr;
		  
		  //   helper indices
		  int *d_diagonalIdx;
		  int *d_kIDinCSRval;
	  } csrMat;

	  int m_node_threadsPerBlock;
	  int m_node_blocksPerGrid;
	  int m_tet_threadsPerBlock;
	  int m_tet_blocksPerGrid;

};