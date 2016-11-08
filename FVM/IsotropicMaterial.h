#pragma once
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen\Dense>
#include <Eigen\Core>
#include <Eigen\svd>
#include <Eigen\Sparse>
#include <vector>
#include <memory>
#include <array>
#include <QTime>

#include "TetMesh.h"
#include "omp.h"

class IsotropicMaterial
{
public:
	IsotropicMaterial();
	virtual ~IsotropicMaterial();

	//virtual float computeEnergy(int tetID, Eigen::Vector3f invariants);
	virtual Eigen::Vector3f computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3f invariants)=0;
	virtual Eigen::Matrix3f computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3f invariants)=0;
	
	// Use separable elastic strain Energy: See [Xu, Sin, Zhu, Barbic 2015] paper Section 4.1
	// "Nonlinear Material Design Using Principal Stretches"
	// Here I do not use high-level eigen data structure for later
	// it's easier to tranfer the codes to CUDA version
	void computeEnergy2FhatGradient(int tetID, const float *Fhats, float *gradient);
	void computeEnergy2FhatHessian(int tetID, const float *Fhats, float *hessian);

	Eigen::MatrixXf computeInnerForcesfromFhats();
	Eigen::MatrixXf computeInnerForcesfromFhats(int num_Threads);
	// compute stiffness Matrix
	Eigen::MatrixXf computeStiffnessMatrix(int tetID);
	Eigen::SparseMatrix<float> computeGlobalStiffnessMatrix();
	Eigen::SparseMatrix<float> computeGlobalStiffnessMatrix(int num_Threads);

	const float m_eps_singularvalue;

protected:
	
	virtual float fEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float gEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float hEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float dfEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float dgEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float dhEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float ddfEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float ddgEnergy(float x, const float & mu, const float & lambda) = 0;
	virtual float ddhEnergy(float x, const float & mu, const float & lambda) = 0;

	void allocateGlobalStiffnessMatrix();
	std::vector<int> m_reserveSize;
	Eigen::SparseMatrix<float> m_globalK;

	std::vector<float> m_mus;
	std::vector<float> m_lambdas;
	std::shared_ptr<TetMesh> m_tetModel;
	Eigen::MatrixXf m_Fhats;
	Eigen::MatrixXf m_Us;
	Eigen::MatrixXf m_Vs;
	Eigen::MatrixXf m_Invariants;

private:
	//see [Teran. 2004], compute F_hat and make sure U,V are real rotation matrix.
	void computeSVD33modified(Eigen::Matrix3f F, Eigen::Vector3f &S, Eigen::Matrix3f &U, Eigen::Matrix3f &V);
	void computeFhatsInvariants();
	void computeFhatsInvariants(int num_Threads);   // openmp-parallel  
	void computeFhats();

	// compute dP/dF using [Teran  2005] method
	Eigen::MatrixXf computeDP2DF(int tetID);
	Eigen::Matrix3f restoreMatrix33fromTeranVector(Eigen::VectorXf v);
	const std::array<int, 9> m_matrix33fromTeran; // transfer teran's order to 3*3 matrix row major

	// compute the diagonal entries of the dP(Fhat)/dFij
	void computeDPFhat2DFij(const float *U, const float *V, const float * hessian, int i, int j, float *dPFhatdFij_diagonal);

	// Different method: see [Xu et al. 2015] Section3.2 equation (5)
	void computeDPDFij(const float *U, const float *Fhat, const float *V, const float *PFhats, const float *hessian, int i, int j, float *dPdFij);
	void computeDP2DF(int tetID, const float *U, const float *Fhat, const float *V, float *dPdF);

	Eigen::Matrix3f helperMatDiagonalMat(Eigen::Matrix3f A, const float *diagonal, Eigen::Matrix3f B);


	//for test
	QTime m_timeTest;
};

