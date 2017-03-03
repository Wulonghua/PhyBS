#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include "CUDA_MATH.h"       // helper matrix operations 
#include "helper_math.h" // helper vector(float3) operations
#include "CUDALinearSolvers.h"
#include "GlobalHelpers.h"


class CUDAInterface
{
public:
	CUDAInterface(int num_nodes, const float *nodes, const float *restPoses, const int *constraintsMask,
		int num_tets, const int *tets, const float youngs, const float nu, const float density, 
		int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval);
	~CUDAInterface();

	void doBackEuler(float *hostNode);
	void computeInnerforces();
	void computeGlobalStiffnessMatrix();
	void updateNodePositions(const float *hostNode);
	void reset();


	void doDescentOpt(float h, int iterations, bool isUpdateH); // h: timestep

private:
	// Descent Optimize 
	void updateDescentOptOneIter(float alpha, float h, bool isUpdateH);
	float sumEnergy();

	  int n_nodes;
	  int n_tets;
	  int *d_tets;
	  int *d_constraintsMask;  
	  float	*d_nodes;
	  float *d_nodes_next;
	  float *d_nodes_last;
	  float *d_nodes_old;
	  float *d_nodes_explicit;

	  float *d_restPoses;
	  float *d_velocities;
	  float *d_velocities_last;
	  float *d_masses;
	  float *d_gravities;
	  float *d_masses_scaled;    // (1+dumping_alpha*timestep)*M
	  float *d_ANs;
	  float *d_vols;			// tet's volumn
	  float *d_Dm_inverses;
	  float *d_mus;
	  float *d_lambdas;
	  float *d_Fhats;
	  float *d_Us;
	  float *d_Vs;
	  float *d_Hd;				// Hessian diagonal
	  float *d_Energy;

	  // first stroe for global stiffness matrix then for LHS of the linear system.
	  // if jacobi iterative method is used, then store for the C matrix Part (see [Wang 2015]) 
	  struct CSRmatrix csrMat;

	  float *d_b;  // first store for forces then for the RHS of the linear system.

	  int *d_kIDinCSRval;

	  int m_node_threadsPerBlock;
	  int m_node_blocksPerGrid;
	  int m_tet_threadsPerBlock;
	  int m_tet_blocksPerGrid;

	  float m_dumpingAlpha;
	  float m_dumpingBelta;
	  float m_timestep;

	  
	  CUDALinearSolvers * cuLinearSolver;

	  // for Descent Optimize Method
	  float		m_alpha;

	  float		m_energy;
	  float		m_energy_old;

	  float		m_rho;
	  float		m_omega;

	  int		m_profile_k[3];
	  float		m_profile_v[3];
	  int		m_profile_n;
};