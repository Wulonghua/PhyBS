#pragma once
#ifdef __INTELLISENSE__   // just to remove the red underline waring (NVIDIA still has not solved the problem)
float atomicAdd(float *address, float val);
#endif


#include "CUDAInterface.h"

__host__ __device__ 
void computeEnergy2FhatGradient(float mu, float lambda, float *Fhats, float *gradient)
{
	float tmp = lambda *  (log(Fhats[0]) + log(Fhats[1]) + log(Fhats[2])) - mu;

	for (int i = 0; i < 3; ++i)
		gradient[i] = mu*Fhats[i] + tmp / Fhats[i];
}

__host__ __device__
void computeEnergy2FhatHessian(float mu, float lambda, const float *Fhats, float *hessian)
{

	float tmp = lambda * (log(Fhats[0]) + log(Fhats[1]) + log(Fhats[2]));
	float tmp1 = mu + lambda - tmp;
	float tmp2 = tmp - mu;

	float inv_lambda1 = 1 / Fhats[0];
	float inv_lambda2 = 1 / Fhats[1];
	float inv_lambda3 = 1 / Fhats[2];

	// hessian(1,1)
	hessian[0] = mu + tmp1 * inv_lambda1 * inv_lambda1;
	// hessian(2,2)
	hessian[4] = mu + tmp1 * inv_lambda2 * inv_lambda2;
	// hessian(3,3)
	hessian[8] = mu + tmp1 * inv_lambda3 * inv_lambda3;

	// hessian(1,2) = hessian(2,1)
	hessian[1] = hessian[3] = (tmp1 + tmp2) * inv_lambda1 * inv_lambda2;

	// hessian(1,3) = hessian(3,1)
	hessian[2] = hessian[6] = (tmp1 + tmp2) * inv_lambda1 * inv_lambda3;

	// hessian(2,3) = hessian(3,2)
	hessian[5] = hessian[7] = (tmp1 + tmp2) * inv_lambda2 * inv_lambda3;
}

__host__ __device__
void computeDPDFij(float *U, float *Fhats, float *V, float *PFhats, float *hessian, int i, int j, float *dPdFij)
{
	float wU[9] = {0.0};
	float wVT[9] = {0.0};

	float A[4], b[2], x[2];
	
	int kl[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
	int k, l;

	for (int kli = 0; kli < 3; ++kli)
	{
		k = kl[kli][0];
		l = kl[kli][1];
		A[0] = A[3] = Fhats[l];
		A[1] = A[2] = Fhats[k];
		b[0] = U[3 * i + k] * V[3 * j + l];
		b[1] = -U[3 * i + l] * V[3 * j + k];

		cuMath::solveLinearSystem22(A, b, x);

		wU[k*3 +l] = x[0];
		wU[l*3 +k] = -x[0];

		wVT[k*3 + l] = x[1];
		wVT[l*3 + k] = -x[1];
	}

	float dUdFij[9], dVTdFij[9];
	cuMath::Matrix_Product_3(U, wU, dUdFij);
	cuMath::Matrix_Product_T_3(wVT, V, dVTdFij);


	float dPFhatdFij[3];
	for (int k = 0; k < 3; ++k)
	{
		dPFhatdFij[k] = hessian[k] * (U[i * 3] * V[j * 3]) +
			hessian[k + 3] * (U[i * 3 + 1] * V[j * 3 + 1]) +
			hessian[k + 6] * (U[i * 3 + 2] * V[j * 3 + 2]);
	}

	float dPdFij1[9], dPdFij2[9], dPdFij3[9];
	cuMath::Matrix_Diagonal_Matrix_T_Product_3(dUdFij, PFhats, V, dPdFij1);
	cuMath::Matrix_Diagonal_Matrix_T_Product_3(U, dPFhatdFij, V, dPdFij2);
	cuMath::Matrix_Diagonal_Matrix_Product_3(U, PFhats, dVTdFij, dPdFij3);

	for (int k = 0; k < 9; ++k)
	{
		dPdFij[k] = dPdFij1[k] + dPdFij2[k] + dPdFij3[k];
	}

}


__host__ __device__ 
void compute_dfdx(float *U, float *Fhat, float *V, float *AN, float *Dm_inv, float mu, float lambda, float *dfdx)
{

	/*********************************Compute dP2dF******************************************/
	float Ftildes[3];
	// perturbation: handle degenerated cases: see the paragraph between equation (9) and equation (10)
	// in paper [Xu et al. 2015]
	bool isPerturbed = true;
	float eps_singularvalue = 1e-6;
	// attention!!! check to make sure Fhats are already sorted in descending order.
	if (Fhat[0] - Fhat[2] < 2 * eps_singularvalue)
	{
		Ftildes[2] = Fhat[2];
		Ftildes[1] = Ftildes[2] + eps_singularvalue;
		Ftildes[0] = Ftildes[1] + eps_singularvalue;
	}
	else // Fhats[0] - Fhats[2] >= 2 * m_eps_singularvalue
	{
		if ((Fhat[0] - Fhat[1] < eps_singularvalue) && (Fhat[1] - Fhat[2] >= eps_singularvalue))
		{
			Ftildes[2] = Fhat[2];
			Ftildes[1] = Fhat[0] - eps_singularvalue;
			Ftildes[0] = Fhat[0];
		}
		else if ((Fhat[0] - Fhat[1] >= eps_singularvalue) && (Fhat[1] - Fhat[2] < eps_singularvalue))
		{
			Ftildes[2] = Fhat[2];
			Ftildes[1] = Fhat[2] + eps_singularvalue;
			Ftildes[0] = Fhat[0];
		}
		else
		{
			Ftildes[2] = Fhat[2];
			Ftildes[1] = Fhat[1];
			Ftildes[0] = Fhat[0];
			isPerturbed = false;
		}
	}

	float Fnew[9], Unew[9], Fhatnew[3], Vnew[9];
	if (isPerturbed)
	{
		cuMath::Matrix_Diagonal_Matrix_T_Product_3(U, Ftildes, V, Fnew);
		cuMath::Matrix_SVD_3(Fnew, Unew, Fhatnew, Vnew);
	}
	else
	{
		memcpy(Unew, U, sizeof(float) * 9);
		memcpy(Vnew, V, sizeof(float) * 9);
		memcpy(Fhatnew, Fhat, sizeof(float) * 3);
	}

	float PFhat[3];
	computeEnergy2FhatGradient(mu, lambda, Fhatnew, PFhat);

	float hessian[9];
	computeEnergy2FhatHessian(mu, lambda, Fhatnew, hessian);

	float dPdF[81];
	float dPdFij[9];
	int Fid;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Fid = 3 * i + j;
			// transpose Unew and Vnew cause they are in stored in column-major matrix (Eigen default setting)
			computeDPDFij(Unew, Fhatnew, Vnew, PFhat, hessian, i, j, dPdFij);
			// copy dPdFij to dPdF
			for (int k = 0; k < 9; ++k)
			{
				dPdF[k * 9 + Fid] = dPdFij[k];
			}
		}
	}
	/********************************************************************/

	float dGdF[81];
	float BT[9];
	cuMath::Matrix_Inverse_3(AN, BT);
	cuMath::Matrix_Product(BT, dPdF, dGdF, 3,3,9);
	cuMath::Matrix_Product(BT, dPdF + 27, dGdF + 27, 3, 3, 9);
	cuMath::Matrix_Product(BT, dPdF + 54, dGdF + 54, 3, 3, 9);

	float dFdx[108] = { 0.0 };
	dFdx[0]  = dFdx[37] = dFdx[74] = -Dm_inv[0] - Dm_inv[3] - Dm_inv[6];
	dFdx[12] = dFdx[49] = dFdx[86] = -Dm_inv[1] - Dm_inv[4] - Dm_inv[7];
	dFdx[24] = dFdx[61] = dFdx[98] = -Dm_inv[2] - Dm_inv[5] - Dm_inv[8];

	dFdx[3]  = dFdx[40] = dFdx[77]  = Dm_inv[0];
	dFdx[15] = dFdx[52] = dFdx[89]  = Dm_inv[1];
	dFdx[27] = dFdx[64] = dFdx[101] = Dm_inv[2];

	dFdx[6]  = dFdx[43] = dFdx[80]  = Dm_inv[3];
	dFdx[18] = dFdx[55] = dFdx[92]  = Dm_inv[4];
	dFdx[30] = dFdx[67] = dFdx[104] = Dm_inv[5];

	dFdx[9]  = dFdx[46] = dFdx[83]  = Dm_inv[6];
	dFdx[21] = dFdx[58] = dFdx[95]  = Dm_inv[7];
	dFdx[33] = dFdx[70] = dFdx[107] = Dm_inv[8];

	float dGdx[108];
	cuMath::Matrix_Product(dGdF, dFdx, dGdx, 9, 9, 12);

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 12; ++j)
		{
			dfdx[i * 12 + j] = -dGdx[36 * i + j] - dGdx[36 * i + j + 12] - dGdx[36 * i + j + 24];
		}
	}

	int convert_idx[9] = { 0, 3, 6, 1, 4, 7, 2, 5, 8 };
	for (int i = 0; i < 9; ++i)
	{
		memcpy(dfdx + (i + 3) * 12, dGdx + convert_idx[i] * 12, sizeof(float) * 12);
	}
}

__global__ void computeGlobalStiffnessMatrix(float *Us, float *Fhats, float *Vs, int *tets, int num_tets, float *ANs, float *Dm_invs,
											 float *mus, float *lambdas, int *kIDinCSRval, float *gK)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;

	int Ki, Kj, gKi, gKj;
	float dfdx[144];
	compute_dfdx(Us + 9 * i, Fhats + 3 * i, Vs + 9 * i, ANs + 9 * i, Dm_invs + 9 * i, mus[i], lambdas[i], dfdx);

	for (int fi = 0; fi < 4; ++fi)
	{
		for (int fj = 0; fj < 3; ++fj)
		{
			Ki = fi * 3 + fj;
			for (int ni = 0; ni < 4; ++ni)
			{
				for (int nj = 0; nj < 3; ++nj)
				{
					Kj = ni * 3 + nj;
					atomicAdd(gK + kIDinCSRval[Ki + ni] + nj, -dfdx[Ki * 12 + Kj]);
				}
			}
		}
	}

	//t_globalK += m_timeTest.restart();
}

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

__global__ void computeInnerForces(float *nodes, int *tets, int num_tets,
								   float *Dm_inv, float *ANs, float *mus, float *lambdas,
								   float *Us, float *Fhats, float *Vs, float *inner_forces)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;
	float Phat[3], Fhat[3], P[9], F[9], U[9], V[9], forces[9];

	/*******************Compute Fhat************************/
	// compute Ds
	float *pos0 = nodes + 3 * tets[i * 4 + 0];
	float *pos1 = nodes + 3 * tets[i * 4 + 1];
	float *pos2 = nodes + 3 * tets[i * 4 + 2];
	float *pos3 = nodes + 3 * tets[i * 4 + 3];

	float Ds[9];
	for (int k = 0; k < 3; ++k)
	{
		Ds[0 + k * 3] = *(pos1 + k) - *(pos0 + k);
		Ds[1 + k * 3] = *(pos2 + k) - *(pos0 + k);
		Ds[2 + k * 3] = *(pos3 + k) - *(pos0 + k);
	}

	// compute deformation gradient
	cuMath::Matrix_Product_3(Ds, Dm_inv + 9 * i, F);

	//svd
	cuMath::Matrix_SVD_3(F, U, Fhat, V);
	memcpy(Us + 9 * i, U, sizeof(float) * 9);
	memcpy(Fhats + 3*i, Fhat, sizeof(float)*3);
	memcpy(Vs+ 9*i, V,sizeof(float)*9);

	computeEnergy2FhatGradient(mus[i], lambdas[i], Fhat, Phat);

	cuMath::Matrix_Diagonal_Matrix_T_Product_3(U, Fhats, V, P);

	cuMath::Matrix_Product_3(P, ANs + 9 * i, forces);

	float *f0 = inner_forces + 3 * tets[i * 4 + 0];
	float *f1 = inner_forces + 3 * tets[i * 4 + 1];
	float *f2 = inner_forces + 3 * tets[i * 4 + 2];
	float *f3 = inner_forces + 3 * tets[i * 4 + 3];

	atomicAdd(f0, -forces[0] - forces[1] - forces[2]);
	atomicAdd(f0 + 1, -forces[3] - forces[4] - forces[5]);
	atomicAdd(f0 + 2, -forces[6] - forces[7] - forces[8]);

	atomicAdd(f1, forces[0]);
	atomicAdd(f1 + 1, forces[3]);
	atomicAdd(f1 + 2, forces[6]);

	atomicAdd(f2, forces[1]);
	atomicAdd(f2 + 1, forces[4]);
	atomicAdd(f2 + 2, forces[7]);

	atomicAdd(f3, forces[2]);
	atomicAdd(f3 + 1, forces[5]);
	atomicAdd(f3 + 2, forces[7]);
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

CUDAInterface::CUDAInterface(int num_nodes, const float *nodes, int num_tets, const int *tets,
	const float youngs, const float nu, const float density,
	int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval):
n_nodes(num_nodes), n_tets(num_tets)
{
	m_node_threadsPerBlock = 64;
	m_node_blocksPerGrid = (num_nodes + m_node_threadsPerBlock - 1) / m_node_threadsPerBlock;
	m_tet_threadsPerBlock = 64;
	m_tet_blocksPerGrid = (num_tets + m_tet_threadsPerBlock - 1) / m_tet_threadsPerBlock;

	csrMat.nnz = csr_nnz;
	csrMat.m = csr_m;

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
	cudaMalloc((void **)&csrMat.d_Valptr, sizeof(float)*csr_nnz);
	cudaMalloc((void **)&csrMat.d_Colptr, sizeof(int)*csr_nnz);
	cudaMalloc((void **)&csrMat.d_Rowptr, sizeof(int)*csr_m);
	cudaMalloc((void **)&csrMat.d_diagonalIdx, sizeof(int)*n_nodes * 3);
	cudaMalloc((void **)&csrMat.d_kIDinCSRval, sizeof(int)*n_tets * 48);

    //Copy data from host to device
	cudaMemcpy(d_tets, tets, sizeof(int) * 4 * n_tets, cudaMemcpyHostToDevice);
	cudaMemcpy(d_nodes, nodes, sizeof(float) * 3 * n_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Valptr, csr_val, sizeof(float)*csr_nnz, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Colptr, csr_col, sizeof(float)*csr_nnz, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Rowptr, csr_row, sizeof(float)*csr_m, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_diagonalIdx, csr_diagonalIdx, sizeof(int)*n_nodes * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_kIDinCSRval, csr_kIDinCSRval, sizeof(int)*n_tets * 48, cudaMemcpyHostToDevice);

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

