#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "CUDA_MATH.h"


class CUDAInterface
{
public:
	CUDAInterface(int num_nodes, const float *nodes, int num_tets, const int *tets, const float youngs, const float nu);
	~CUDAInterface();


private:
	  int n_nodes;
	  int n_tets;
	  int *d_tets;
	  float	*d_nodes;
	  float *d_ANs;
	  float *d_Dm_inverses;
	  float *d_mus;
	  float *d_lambdas;
	  float *d_Fhats;
	  float *d_Us;
	  float *d_Vs;

	  int m_node_threadsPerBlock;
	  int m_node_blocksPerGrid;
	  int m_tet_threadsPerBlock;
	  int m_tet_blocksPerGrid;

};