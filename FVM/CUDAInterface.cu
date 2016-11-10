#pragma once
#ifdef __INTELLISENSE__   // just to remove the red underline waring (NVIDIA still has not solved the problem)
float atomicAdd(float *address, float val);
#endif


#include "CUDAInterface.h"

__global__ void compute_DmInv_ANs_Mass(float *nodes, int *tets, int num_tets, float density, float *Dm_inv, float *masses, float *ANs)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;
	
	// compute Dm_Inverses
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

	// compute node' mass
	float m = fabsf(cuMath::Matrix_Determinant_3(Dm)) / 24.0 * density;
	for (int k = 0; k < 4; ++k)
	{
		atomicAdd(masses+tets[i * 4 + k], m);
	}


	// compute ANs
	int edge[3][4] = { { 1, 0, 2, 3 }, { 2, 0, 3, 1 }, { 3, 0, 1, 2 } };
	float3 p[4], v1, v2, v3, vtmp;
	float ANs_[9];
	p[0] = make_float3(*pos0, *(pos0 + 1), *(pos0 + 2));
	p[1] = make_float3(*pos1, *(pos1 + 1), *(pos1 + 2));
	p[2] = make_float3(*pos2, *(pos2 + 1), *(pos2 + 2));
	p[3] = make_float3(*pos3, *(pos3 + 1), *(pos3 + 2));

	for (int k = 0; k < 3; ++k)
	{
		v1 = p[edge[k][1]] - p[edge[k][0]];
		v2 = p[edge[k][2]] - p[edge[k][0]];
		v3 = p[edge[k][3]] - p[edge[k][0]];
		vtmp = -(cross(v1, v2) + cross(v2, v3) + cross(v3, v1)) / 6.0;
		ANs_[k] = vtmp.x;
		ANs_[k + 3] = vtmp.y;
		ANs_[k + 6] = vtmp.z;
	}

	memcpy(ANs + 9 * i, ANs_, sizeof(float) * 9);


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
	//printf("v: %f\n",devicePtr[7199]);
}

CUDAInterface::CUDAInterface(int num_nodes, const float *nodes, int num_tets, const int *tets, const float youngs, const float nu, const float density) :
n_nodes(num_nodes), n_tets(num_tets)
{
	m_node_threadsPerBlock = 64;
	m_node_blocksPerGrid = (num_nodes + m_node_threadsPerBlock - 1) / m_node_threadsPerBlock;
	m_tet_threadsPerBlock = 64;
	m_tet_blocksPerGrid = (num_tets + m_tet_threadsPerBlock - 1) / m_tet_threadsPerBlock;

	// allocate device memory
	cudaMalloc((void **)&d_tets, sizeof(int) * 4 * n_tets);
	cudaMalloc((void **)&d_nodes, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_masses, sizeof(float) * n_nodes);
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

	//initialize weights to 0
	cudaMemset(d_masses, 0, sizeof(float)*n_nodes);

	// currently using homogenous material
	float mu = 0.5 * youngs / (1 + nu);
	float lambda = mu * 2 / (1 - 2 * nu);

	//initialize d_mus and d_lambda
	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_mus, mu, n_tets);
	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_lambdas, lambda, n_tets);
	
	// compute d_Dm_inverses, d_masses, and d_ANs
	compute_DmInv_ANs_Mass << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_nodes, d_tets,
																 n_tets, density, d_Dm_inverses,
																				d_masses, d_ANs);	
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

