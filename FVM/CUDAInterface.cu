#pragma once

#include "CUDAInterface.h"

__global__ void computeDmInverse(float *nodes, int *tets, int num_tets, float *Dm_inv)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;
	
	float *pos0 = nodes + 3 * tets[i * 4 + 0];
	float *pos1 = nodes + 3 * tets[i * 4 + 1];
	float *pos2 = nodes + 3 * tets[i * 4 + 2];
	float *pos3 = nodes + 3 * tets[i * 4 + 3];

	float Dm[9],Dm_[9];
	for (int k = 0; k < 3; ++k)
	{
		Dm[0+k*3] = *(pos1+k)-*(pos0+k);
		Dm[1+k*3] = *(pos2+k)-*(pos0+k);
		Dm[2+k*3] = *(pos3+k)-*(pos0+k);
	}
	cuMath::Matrix_Inverse_3(Dm, Dm_);
	memcpy(Dm_inv + 9 * i, Dm_, sizeof(float) * 9);
}

__global__ void setDeviceArray(float *devicePtr, float v, int num)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num)	return;

	devicePtr[i] = v;
}


// this is helper function that is used for NSIGHT debug.
__global__ void checkData(float *devicePtr)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
}

CUDAInterface::CUDAInterface(int num_nodes, const float *nodes, int num_tets, const int *tets, const float youngs, const float nu) :
n_nodes(num_nodes), n_tets(num_tets)
{
	m_node_threadsPerBlock = 64;
	m_node_blocksPerGrid = (num_nodes + m_node_threadsPerBlock - 1) / m_node_threadsPerBlock;
	m_tet_threadsPerBlock = 64;
	m_tet_blocksPerGrid = (num_tets + m_tet_threadsPerBlock - 1) / m_tet_threadsPerBlock;

	// allocate device memory
	cudaMalloc((void **)&d_tets, sizeof(int) * 4 * n_tets);
	cudaMalloc((void **)&d_nodes, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_ANs, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_Dm_inverses, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_mus, sizeof(float) * n_tets);
	cudaMalloc((void **)&d_lambdas, sizeof(float) * n_tets);
	cudaMalloc((void **)&d_Fhats, sizeof(float) * 3 * n_tets);
	cudaMalloc((void **)&d_Us, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_Vs, sizeof(float) * 9 * n_tets);

    //Copy data from host to device
	cudaMemcpy(d_tets, tets, sizeof(int) * 4 * n_tets, cudaMemcpyHostToDevice);
	cudaMemcpy(d_nodes, nodes, sizeof(float) * 3 * n_nodes, cudaMemcpyHostToDevice);

	// currently using homogenous material
	float mu = 0.5 * youngs / (1 + nu);
	float lambda = mu * 2 / (1 - 2 * nu);

	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_mus, mu, n_tets);
	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_lambdas, lambda, n_tets);
	
}

CUDAInterface::~CUDAInterface()
{
	cudaFree(d_tets);
	cudaFree(d_nodes);
	cudaFree(d_ANs);
	cudaFree(d_Dm_inverses);
	cudaFree(d_mus);
	cudaFree(d_lambdas);
	cudaFree(d_Fhats);
	cudaFree(d_Us);
	cudaFree(d_Vs);
}

