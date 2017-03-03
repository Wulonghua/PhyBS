#pragma once
#ifdef __INTELLISENSE__   // just to remove the red underline waring (NVIDIA has not solved the problem)
float atomicAdd(float *address, float val);
#endif


#include "CUDAInterface.h"


__host__ __device__ 
void computeEnergy2FhatGradient(float mu, float lambda, float *Fhats, float *gradient)
{
	float tmp = lambda *  (std::log(Fhats[0]) + std::log(Fhats[1]) + std::log(Fhats[2])) - mu;

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
	cuMath::Matrix_Transose_3(AN, BT);
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

__host__ __device__ void computeFirstPiolaStress(float *F, float mu, float lambda, float *P)
{
	float F_invT[9];
	cuMath::Matrix_Inverse_Transpose_3(F, F_invT);

	float F_determinant = cuMath::Matrix_Determinant_3(F);

	for (int i = 0; i < 9; ++i)
	{
		P[i] = (F[i] - F_invT[i]) * mu + lambda * log(F_determinant) * F_invT[i];
	}
}

__global__ void computeGlobalK(float *Us, float *Fhats, float *Vs, int *tets, int num_tets, float *ANs, float *Dm_invs, float *mus, float *lambdas, int *kIDinCSRval, float *gK)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;

	int Ki, Kj,KIdx;
	int tetInitID = i * 48;
	float dfdx[144];
	compute_dfdx(Us + 9 * i, Fhats + 3 * i, Vs + 9 * i, ANs + 9 * i, Dm_invs + 9 * i, mus[i], lambdas[i], dfdx);

	//for (int ti = 0; ti < 12; ++ti)
	//{
	//	for (int tj = 0; tj < 12; ++tj)
	//	{
	//		printf(" %f", dfdx[ti * 12 + tj]);
	//	}
	//	printf("\n");
	//}

	for (int fi = 0; fi < 4; ++fi)
	{
		for (int fj = 0; fj < 3; ++fj)
		{
			Ki = fi * 3 + fj;
			for (int ni = 0; ni < 4; ++ni)
			{
				KIdx = Ki * 4 + ni;
				for (int nj = 0; nj < 3; ++nj)
				{
					Kj = ni * 3 + nj;
					atomicAdd(gK + kIDinCSRval[tetInitID+KIdx] + nj, -dfdx[Ki * 12 + Kj]);
				}
			}
		}
	}

	//t_globalK += m_timeTest.restart();
}

__global__ void compute_DmInv_ANs_Mass(float *nodes, int *tets, int num_tets, float density, float *Dm_inv, float *masses, float *ANs, float *vols)
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

	// compute tet's volumn and node' mass
	float v = fabsf(cuMath::Matrix_Determinant_3(Dm)) / 6.0;
	vols[i] = v;
	float m = v * 0.25 * density;
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

__global__ void addFixContraintsAndGravity(float *nodes, float *rest, float *gravity, float *inner_forces, int *mask, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	if (mask[i] > 0)
	{
		inner_forces[3 * i + 0] = 1e8 *(rest[3 * i + 0] - nodes[3 * i + 0]);
		inner_forces[3 * i + 1] = 1e8 *(rest[3 * i + 1] - nodes[3 * i + 1]);
		inner_forces[3 * i + 2] = 1e8 *(rest[3 * i + 2] - nodes[3 * i + 2]);
		
	}
	else
	{
		inner_forces[3 * i + 1] += gravity[i];
	}

}

__global__ void computeElasticForces(float *nodes, int *tets, int num_tets,
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

	//memset(F, 0, sizeof(float) * 9);
	//F[0] = F[4] = F[8] = 1.0;

	//for (int i = 0; i < 3; ++i)
	//{
	//	for (int j = 0; j < 3; ++j)
	//	{
	//		printf(" %.16f", F[i * 3 + j]);
	//	}
	//	printf("\n");
	//}

	//svd
	cuMath::Matrix_SVD_3(F, U, Fhat, V);

	memcpy(Us + 9 * i, U, sizeof(float) * 9);
	memcpy(Fhats + 3*i, Fhat, sizeof(float)*3);
	memcpy(Vs+ 9*i, V,sizeof(float)*9);

	computeEnergy2FhatGradient(mus[i], lambdas[i], Fhat, Phat);

	cuMath::Matrix_Diagonal_Matrix_T_Product_3(U, Phat, V, P);

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
	atomicAdd(f3 + 2, forces[8]);
}

__global__ void setRHSofLinearSystem(float *m, float *v, float t, int n, float *b)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	//if (i == 0)
	//{
	//	printf("%.10f %.10f\n", b[3 * i + 0], v[3*i+0]);
	//	printf("%.10f %.10f\n", b[3 * i + 1], v[3*i+1]);
	//	printf("%.10f %.10f\n", b[3 * i + 2], v[3*i+2]);
	//}

	b[3 * i + 0] = b[3 * i + 0] * t + m[i] * v[3 * i + 0];
	b[3 * i + 1] = b[3 * i + 1] * t + m[i] * v[3 * i + 1];
	b[3 * i + 2] = b[3 * i + 2] * t + m[i] * v[3 * i + 2];

	//if (i == 0)
	//{
	//	printf("%f\n", b[3 * i + 0]);
	//	printf("%f\n", b[3 * i + 1]);
	//	printf("%f\n", b[3 * i + 2]);
	//}

}

// do not use this function which has unsolved strange bug
__global__ void setLHSofLinearSystem(float *m, float t, float dumpingAlpha, int n, float *Valptr, int *diagonalIdx)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;
	float tmp = m[i] * (1 + dumpingAlpha *t);
	//for (int j = 0; j < 3; ++j)
	//{
	//	id = i * 3 + j;
	//	Valptr[diagonalIdx[id]] += m[i] * (1 + dumpingAlpha * t);
	//}

	/*********************Strange bug here.******************************/
	
	if (i == 0)
	{
		Valptr[diagonalIdx[i * 3 + 0]] += tmp;
		Valptr[diagonalIdx[i * 3 + 1]] += tmp;
		Valptr[diagonalIdx[i * 3 + 2]] += tmp;
	}
}

__global__ void setLHSofLinearSystem(float *m_s, int n, float *Valptr, int *diagonalIdx)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	Valptr[diagonalIdx[i]] += m_s[i];
}

__global__ void setMassScaled(float *mass, float * mass_s, float s,int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	mass_s[3 * i + 0] = mass_s[3 * i + 1] = mass_s[3 * i + 2] = mass[i] * s;
}

__global__ void initDescentOpt(float *nodes, float *nodes_explicit, float *v, float *v_last, float h, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n) return;
	int id;
	for (int k = 0; k < 3; ++k)
	{
		id = i * 3 + k;
		nodes_explicit[id] = nodes[id] + v[id] * h;
		nodes[id] = nodes_explicit[id] + 0.2 * (v[id] - v_last[id]) * h;
	}
}

/****************** Descend Optimize Method***********************************************/
__global__ void computeForceKd(bool computeHd, float *nodes, int *tets, int num_tets, float *Dm_invs, float *ANs, float *mus, float *lambdas, float *vols, float *Kd, float *inner_forces,float *Energys)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_tets)	return;

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
	float F[9], P[9], forces[9];
	cuMath::Matrix_Product_3(Ds, Dm_invs + 9 * i, F);
	computeFirstPiolaStress(F, mus[i], lambdas[i], P);

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
	atomicAdd(f3 + 2, forces[8]);

	// compute elastic energy
	float energy, logJ;
	float I1 = 0.0;
	for (int k = 0; k < 9; ++k)
		I1 += F[k] * F[k];
	logJ = log(cuMath::Matrix_Determinant_3(F));
	energy = (0.5 * mus[i] * (I1 - 3) - mus[i] * logJ + 0.5 * lambdas[i] * logJ * logJ)*vols[i];
	energy = MAX(energy, 0.0);
	atomicAdd(Energys+tets[i*4], energy);

	if (!computeHd) return;

	float Fhat[9], U[9], V[9];
	cuMath::Matrix_SVD_3(F, U, Fhat, V);

	float dfdx[144];
	compute_dfdx(U, Fhat, V, ANs + 9 * i, Dm_invs + 9 * i, mus[i], lambdas[i], dfdx);

	for (int n = 0; n < 4; ++n)
	{
		for (int m = 0; m < 3; ++m)
		{
			int k = n * 3 + m;
			atomicAdd(Kd + tets[i * 4 + n] * 3 + m, -dfdx[k * 12 + k]);
		}
	}
}

// compute delta_q (see CPU version codes for reference) and update nodes_next, also add energy
// b: inner force
__global__ void updateNodesNextIter(float *Hd, float *b , float *masses, float *nodes, float *nodes_explicit, float *v, float h, float *gravitys, int *fixed_mask, float alpha, int num_nodes, float *Energys, float *nodes_next)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_nodes)	return;
	if (fixed_mask[i] > 0) return;

	// adds energys
	int id;
	float m_h2_inv = masses[i] / (h*h);

	for (int k = 0; k < 3; ++k)
	{
		id = i * 3 + k;
		Energys[i] += 0.5 * (nodes[id] - nodes_explicit[id]) * (nodes[id] - nodes_explicit[id]) * m_h2_inv; 
	}
	Energys[i] += 9.8 * nodes[i * 3 + 1] * masses[i];

	float g[3];
	g[0] = (nodes[i * 3    ] - nodes_explicit[i * 3    ]) * m_h2_inv - b[i * 3    ];
	g[1] = (nodes[i * 3 + 1] - nodes_explicit[i * 3 + 1]) * m_h2_inv - b[i * 3 + 1] - gravitys[i];
	g[2] = (nodes[i * 3 + 2] - nodes_explicit[i * 3 + 2]) * m_h2_inv - b[i * 3 + 2];

	for (int k = 0; k < 3; ++k)
	{
		id = i * 3 + k;
		nodes_next[id] = nodes[id] - g[k] / (Hd[id]+m_h2_inv);
	}
}

__global__ void updateChebyshev(float omega, int* fixed_mask, int num_nodes, float *nodes_next, float *nodes_last)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_nodes)	return;
	if (fixed_mask[i] > 0) return;

	int id;
	for (int k = 0; k < 3; ++k)
	{
		id = 3 * i + k;
		nodes_next[id] = omega * (nodes_next[id] - nodes_last[id]) + nodes_last[id];
	}
}

__global__ void updateVelocity(int num_nodes, float inv_h, int *fixed_mask, float *nodes, float *nodes_explicit, float *v, float *v_last)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num_nodes)	return;
	if (fixed_mask[i] > 0) return;

	int id;
	for (int k = 0; k < 3; ++k)
	{
		id = 3 * i + k;
		v_last[id] = v[id];
		v[id] += (nodes[id] - nodes_explicit[id]) * inv_h;
	}
}
/*******************************************************************************************/
CUDAInterface::CUDAInterface(int num_nodes, const float *nodes, const float *restPoses, const int *constraintsMask,
	int num_tets, const int *tets, const float youngs, const float nu, const float density,
	int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval):
	n_nodes(num_nodes), n_tets(num_tets), m_dumpingAlpha(0.02), m_dumpingBelta(0.02), m_timestep(0.03)
{
	m_node_threadsPerBlock = 64;
	m_node_blocksPerGrid = (num_nodes + m_node_threadsPerBlock - 1) / m_node_threadsPerBlock;
	m_tet_threadsPerBlock = 64;
	m_tet_blocksPerGrid = (num_tets + m_tet_threadsPerBlock - 1) / m_tet_threadsPerBlock;

	csrMat.nnz = csr_nnz;
	csrMat.m = csr_m;    // matrix's row number

	cuLinearSolver = new CUDALinearSolvers(csr_m, csr_nnz);

	// allocate device memory
	cudaMalloc((void **)&d_tets, sizeof(int) * 4 * n_tets);
	cudaMalloc((void **)&d_nodes, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_nodes_next, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_nodes_last, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_nodes_old, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_nodes_explicit, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_restPoses, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_velocities, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_velocities_last, sizeof(float) * 3 * n_nodes);
	cudaMalloc((void **)&d_masses, sizeof(float) * n_nodes);
	cudaMalloc((void **)&d_gravities, sizeof(float) * n_nodes);
	cudaMalloc((void **)&d_masses_scaled, sizeof(float) * n_nodes*3);
	cudaMalloc((void **)&d_constraintsMask, sizeof(int) * n_nodes);
	cudaMalloc((void **)&d_ANs, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_Dm_inverses, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_mus, sizeof(float) * n_tets);
	cudaMalloc((void **)&d_lambdas, sizeof(float) * n_tets);
	cudaMalloc((void **)&d_vols, sizeof(float) * n_tets);
	cudaMalloc((void **)&d_Fhats, sizeof(float) * 3 * n_tets);
	cudaMalloc((void **)&d_Us, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&d_Vs, sizeof(float) * 9 * n_tets);
	cudaMalloc((void **)&csrMat.d_Valptr, sizeof(float)*csr_nnz);
	cudaMalloc((void **)&csrMat.d_Colptr, sizeof(int)*csr_nnz);
	cudaMalloc((void **)&csrMat.d_Rowptr, sizeof(int)*(csr_m+1));
	cudaMalloc((void **)&csrMat.d_diagonalIdx, sizeof(int)*n_nodes * 3);
	cudaMalloc((void **)&d_kIDinCSRval, sizeof(int)*n_tets * 48);
	cudaMalloc((void **)&d_b, sizeof(float) * n_nodes * 3);
	cudaMalloc((void **)&d_Hd, sizeof(float) * n_nodes * 3);
	cudaMalloc((void **)&d_Energy, sizeof(float) * n_nodes);

    //Copy data from host to device
	cudaMemcpy(d_tets, tets, sizeof(int) * 4 * n_tets, cudaMemcpyHostToDevice);
	cudaMemcpy(d_nodes, nodes, sizeof(float) * 3 * n_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_restPoses, restPoses, sizeof(float) * 3 * n_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(d_constraintsMask, constraintsMask, sizeof(int) * n_nodes, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Valptr, csr_val, sizeof(float)*csr_nnz, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Colptr, csr_col, sizeof(int)*csr_nnz, cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_Rowptr, csr_row, sizeof(int)*(csr_m+1), cudaMemcpyHostToDevice);
	cudaMemcpy(csrMat.d_diagonalIdx, csr_diagonalIdx, sizeof(int)*n_nodes * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(d_kIDinCSRval, csr_kIDinCSRval, sizeof(int)*n_tets * 48, cudaMemcpyHostToDevice);

	//initialize weights to 0
	cudaMemset(d_masses, 0, sizeof(float)*n_nodes);
	cudaMemset(d_velocities, 0, sizeof(float)*n_nodes*3);

	// currently using homogenous material
	float mu = 0.5 * youngs / (1 + nu);
	float lambda = mu * 2 / (1 - 2 * nu);

	//initialize d_mus and d_lambda
	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_mus, mu, n_tets);
	setDeviceArray << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_lambdas, lambda, n_tets);
	cudaDeviceSynchronize();
	// compute d_Dm_inverses, d_masses, and d_ANs
	compute_DmInv_ANs_Mass << < m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_restPoses, d_tets,
																 n_tets, density, d_Dm_inverses,
																				d_masses, d_ANs, d_vols);
	cudaDeviceSynchronize();
	scaleDeviceArray << < m_node_blocksPerGrid, m_node_threadsPerBlock >> >(d_masses, d_gravities, -9.8, n_nodes);
	cudaDeviceSynchronize();
	// set fixed nodes to relative large masses
	setDeviceArray << < m_node_blocksPerGrid, m_node_threadsPerBlock >> > (d_masses, d_constraintsMask, 1e8, n_nodes);
	cudaDeviceSynchronize();
	//std::cout << "mass: " << std::endl;
	//printData << <1, 1 >> >(d_masses, 4, 1);
	//cudaDeviceSynchronize();

	setMassScaled << < m_node_blocksPerGrid, m_node_threadsPerBlock >> >(d_masses, d_masses_scaled, (1+m_dumpingAlpha*m_timestep), n_nodes);
	cudaDeviceSynchronize();

	//std::cout << "scaled mass: " << std::endl;
	//printData << <1, 1 >> >(d_masses_scaled, 4, 3);
	//cudaDeviceSynchronize();

	m_profile_k[0] = 1;
	m_profile_k[1] = 7;
	m_profile_k[2] = 10;
	m_profile_v[0] = 0.96;
	m_profile_v[1] = 0.991;
	m_profile_v[2] = 0.999;
	m_profile_n = 3;
	reset();
}

CUDAInterface::~CUDAInterface()
{
	cudaFree(d_tets);
	cudaFree(d_nodes);
	cudaFree(d_nodes_next);
	cudaFree(d_nodes_last);
	cudaFree(d_nodes_old);
	cudaFree(d_nodes_explicit);
	cudaFree(d_ANs);
	cudaFree(d_Dm_inverses);
	cudaFree(d_mus);
	cudaFree(d_lambdas);
	cudaFree(d_vols);
	cudaFree(d_Fhats);
	cudaFree(d_Us);
	cudaFree(d_Vs);
	cudaFree(d_constraintsMask);
	cudaFree(d_restPoses);
	cudaFree(d_velocities);
	cudaFree(d_velocities_last);
	cudaFree(d_masses);
	cudaFree(d_gravities);
	cudaFree(d_masses_scaled);
	cudaFree(d_b);
	cudaFree(d_Hd);
	cudaFree(d_Energy);
}

void CUDAInterface::updateNodePositions(const float *hostnodes)
{
	cudaMemcpy(d_nodes, hostnodes, sizeof(float) * 3 * n_nodes, cudaMemcpyHostToDevice);
}

void CUDAInterface::reset()
{
	cudaMemcpy(d_nodes, d_restPoses, sizeof(float) * 3 * n_nodes, cudaMemcpyDeviceToDevice);
	cudaMemcpy(d_nodes_last, d_restPoses, sizeof(float) * 3 * n_nodes, cudaMemcpyDeviceToDevice);
	cudaMemset(d_velocities, 0, sizeof(float) * 3 * n_nodes);
	cudaMemset(d_velocities_last, 0, sizeof(float) * 3 * n_nodes);
}

void CUDAInterface::computeInnerforces()
{
	cudaMemset(d_b, 0, sizeof(float)*n_nodes * 3);
	computeElasticForces << <m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_nodes, d_tets, n_tets, d_Dm_inverses, d_ANs, d_mus, d_lambdas, d_Us, d_Fhats, d_Vs, d_b);
	cudaDeviceSynchronize();

	addFixContraintsAndGravity << < m_node_blocksPerGrid, m_node_threadsPerBlock >> >(d_nodes, d_restPoses, d_gravities, d_b, d_constraintsMask, n_nodes);
	cudaDeviceSynchronize();
	//std::cout << "dUs:" << std::endl;
	//printData<< <1, 1 >> >(d_Us, 6, 3);
	//cudaDeviceSynchronize();

	//std::cout << "dVs:" << std::endl;
	//printData << <1, 1 >> >(d_Vs, 6, 3);
	//cudaDeviceSynchronize();

	//std::cout << "Fhats: " << std::endl;
	//printData << <1, 1 >> >(d_Fhats, 6, 1);
	//cudaDeviceSynchronize();
	//printData << <1, 1 >> >(d_b, n_nodes, 3);
	//cudaDeviceSynchronize();

	//printData << < 1, 1 >> >(d_gravities, n_nodes, 1);
	//cudaDeviceSynchronize();
}

void CUDAInterface::computeGlobalStiffnessMatrix()
{
	cudaMemset(csrMat.d_Valptr, 0, sizeof(float)*csrMat.nnz);
	computeGlobalK << <m_tet_blocksPerGrid, m_tet_threadsPerBlock >> >(d_Us, d_Fhats, d_Vs, d_tets, n_tets, d_ANs, d_Dm_inverses, d_mus, d_lambdas, d_kIDinCSRval, csrMat.d_Valptr);
	cudaDeviceSynchronize();
	//printData << <1, 1 >> >(csrMat.d_Valptr, 12, 12);
	//cudaDeviceSynchronize();
}

void CUDAInterface::doBackEuler(float *hostNode)
{

	//std::cout << "velocities: " << std::endl;
	//printData << <1, 1 >> >(d_velocities, 12, 1);
	//cudaDeviceSynchronize();

	//std::cout << "masses: " << std::endl;
	//printData << <1, 1 >> >(d_masses, 1, 5);
	//cudaDeviceSynchronize();
	//std::cout << "nnz:" << std::endl;
	//std::cout << csrMat.nnz << std::endl;

	//std::cout << "idx: " << std::endl;
	//printData << <1, 1 >> >(d_kIDinCSRval, 1, csrMat.nnz);
	//cudaDeviceSynchronize();

	//std::cout << "K: " << std::endl;
	//printData << <1, 1 >> >(csrMat.d_Valptr, n_nodes*3, n_nodes*3);
	//cudaDeviceSynchronize();

	//std::cout << "forces: " << std::endl;
	//printData << <1, 1 >> >(d_b, n_nodes, 3);
	//cudaDeviceSynchronize();

	setRHSofLinearSystem << < m_node_blocksPerGrid, m_node_threadsPerBlock >> >(d_masses, d_velocities, m_timestep, n_nodes, d_b);
	cudaDeviceSynchronize();
	int threadsPerBlock = 64;
	int blocksPerGrid = (csrMat.nnz + threadsPerBlock - 1) / threadsPerBlock;
	scaleDeviceArray<<<blocksPerGrid, threadsPerBlock>>>(csrMat.d_Valptr, m_timestep + m_dumpingBelta, m_timestep, csrMat.nnz);
	cudaDeviceSynchronize();

	setLHSofLinearSystem << <m_node_blocksPerGrid, m_node_threadsPerBlock * 3 >> >(d_masses_scaled, csrMat.m, csrMat.d_Valptr, csrMat.d_diagonalIdx);
	cudaDeviceSynchronize();

	//std::cout << "RHS: " << std::endl;
	//printData << <1, 1 >> >(d_b, n_nodes*3, 1);
	//cudaDeviceSynchronize();

	//std::cout << "LHS: " << std::endl;
	//printData << <1, 1 >> >(csrMat.d_Valptr, n_nodes*3, n_nodes*3);
	//cudaDeviceSynchronize();

	// has bug
	//cuLinearSolver->chebyshevSemiIteration(csrMat.d_Valptr, csrMat.d_Rowptr, csrMat.d_Colptr, csrMat.d_diagonalIdx,d_b, 0.999, &d_velocities);

	//cuLinearSolver->jacobiIteration(csrMat.d_Valptr, csrMat.d_Rowptr, csrMat.d_Colptr, csrMat.d_diagonalIdx, d_b, &d_velocities);

	//cuLinearSolver->conjugateGradient(csrMat.d_Valptr, csrMat.d_Rowptr, csrMat.d_Colptr, d_b, d_velocities);
	cuLinearSolver->directCholcuSolver(csrMat.d_Valptr, csrMat.d_Rowptr, csrMat.d_Colptr, d_b, d_velocities);
	//cudaDeviceSynchronize();
	//std::cout << "velocities: " << std::endl;
	//printData << < 1, 1 >> > (d_velocities,n_nodes*3,1);
	//cudaDeviceSynchronize();
	
	// update positions. 
	cublasSaxpy(cuLinearSolver->getcuBlasHandle(), n_nodes * 3, &m_timestep, d_velocities, 1, d_nodes, 1);
	cudaDeviceSynchronize();

	//std::cout << "p: " << std::endl;
	//printData << < 1, 1 >> > (d_nodes, n_nodes * 3, 1);
	//cudaDeviceSynchronize();

	cudaMemcpy(hostNode, d_nodes, sizeof(float)*n_nodes * 3, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

float CUDAInterface::sumEnergy()
{
	thrust::device_ptr<float> d_e(d_Energy);
	return thrust::reduce(d_e, d_e + n_nodes);
}

void CUDAInterface::updateDescentOptOneIter(float alpha, float h, bool isUpdateH)
{
	if (isUpdateH)
		cudaMemset(d_Hd,	 0, sizeof(float) * n_nodes * 3);
	cudaMemset(d_b,		 0, sizeof(float) * n_nodes); // d_b is force
	cudaMemset(d_Energy, 0, sizeof(float) * n_nodes);

	computeForceKd<<<m_tet_blocksPerGrid, m_tet_threadsPerBlock>>>(isUpdateH, d_nodes, d_tets, n_tets, d_Dm_inverses, d_ANs, d_mus, d_lambdas, d_vols, d_Hd, d_b, d_Energy);
	updateNodesNextIter<<<m_node_blocksPerGrid, m_node_threadsPerBlock>>>(d_Hd, d_b, d_masses, d_nodes, d_nodes_explicit, d_velocities, m_timestep,
		d_gravities, d_constraintsMask, alpha, n_nodes, d_Energy, d_nodes_next);
}

void CUDAInterface::doDescentOpt(float h, int iterations, float *hostNode)
{
	m_alpha = 0.1;
	initDescentOpt<<<m_node_blocksPerGrid, m_node_threadsPerBlock>>>(d_nodes, d_nodes_explicit, d_velocities, d_velocities_last, m_timestep, n_nodes);
	
	for (int k = 0; k < iterations; ++k)
	{
		if (m_alpha < 0.01) break;
		bool isUpdateH = false;
		if (k % 32 == 0)
			isUpdateH = true;

		updateDescentOptOneIter(m_alpha, m_timestep, isUpdateH);

		if (k % 8 == 0)
		{
			m_energy = sumEnergy();

			if (k != 0 && m_energy > m_energy_old)
			{
				m_alpha *= 0.7;
				iterations -= k - 8;
				k -= 1;
				cudaMemcpy(d_nodes_old, d_nodes, sizeof(float)*n_nodes, cudaMemcpyDeviceToDevice);
				continue;
			}
			else
			{
				m_energy_old = m_energy;
				cudaMemcpy(d_nodes_old, d_nodes, sizeof(float)*n_nodes, cudaMemcpyDeviceToDevice);
			}
		}

		if (k == 0)
		{
			m_rho = 0;
			m_omega = 1;
		}
		m_omega = 4 / (4 - m_rho*m_rho*m_omega);
		for (int i = 0; i<m_profile_n; i++)
		{
			if (k == m_profile_k[i] - 1)
			{
				m_rho = 0;
				m_omega = 1;
			}
			if (k == m_profile_k[i])
			{
				m_rho = m_profile_v[i];
				m_omega = 2 / (2 - m_rho*m_rho);
				break;
			}
		}

		if (m_omega != 1)
			updateChebyshev<<<m_node_blocksPerGrid, m_node_threadsPerBlock >>>(m_omega, d_constraintsMask, n_nodes, d_nodes_next, d_nodes_last);

		swap(&d_nodes_last, &d_nodes);
		swap(&d_nodes, &d_nodes_next);
	}

	updateVelocity<<<m_node_blocksPerGrid, m_node_threadsPerBlock>>>(n_nodes, 1.0 / m_timestep, d_constraintsMask, d_nodes, d_nodes_explicit, d_velocities, d_velocities_last);
	cudaMemcpy(hostNode, d_nodes, sizeof(float)*n_nodes * 3, cudaMemcpyDeviceToHost);
}