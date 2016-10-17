#pragma once
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen\Dense>
#include <Eigen\Core>
#include <Eigen\svd>
#include <Eigen\Sparse>
#include <vector>
#include <memory>
#include <array>

#include "TetMesh.h"

class IsotropicMaterial
{
public:
	IsotropicMaterial();
	virtual ~IsotropicMaterial();

	//virtual double computeEnergy(int tetID, Eigen::Vector3d invariants);
	virtual Eigen::Vector3d computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants)=0;
	virtual Eigen::Matrix3d computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants)=0;
	
	// Use separable elastic strain Energy: See [Xu, Sin, Zhu, Barbic 2015] paper Section 4.1
	// "Nonlinear Material Design Using Principal Stretches"
	// Here I do not use high-level eigen data structure for later
	// it's easier to tranfer the codes to CUDA version
	void computeEnergy2FhatGradient(int tetID, const double *Fhats, double *gradient);
	void computeEnergy2FhatHessian(int tetID, const double *Fhats, double *hessian);

	Eigen::MatrixXd computeInnerForcesfromFhats();
	// compute stiffness Matrix
	Eigen::MatrixXd computeStiffnessMatrix(int tetID);
	Eigen::SparseMatrix<double> computeGlobalStiffnessMatrix();

	const double m_eps_singularvalue;

protected:
	
	virtual double fEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double gEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double hEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double dfEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double dgEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double dhEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double ddfEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double ddgEnergy(double x, const double & mu, const double & lambda) = 0;
	virtual double ddhEnergy(double x, const double & mu, const double & lambda) = 0;

	std::vector<double> m_mus;
	std::vector<double> m_lambdas;
	std::shared_ptr<TetMesh> m_tetModel;
	Eigen::MatrixXd m_Fhats;
	Eigen::MatrixXd m_Us;
	Eigen::MatrixXd m_Vs;
	Eigen::MatrixXd m_Invariants;

private:
	//see [Teran. 2004], compute F_hat and make sure U,V are real rotation matrix.
	void computeSVD33modified(Eigen::Matrix3d F, Eigen::Vector3d &S, Eigen::Matrix3d &U, Eigen::Matrix3d &V);
	void computeFhatsInvariants();
	void compouteFhats();

	// compute dP/dF using [Teran  2005] method
	Eigen::MatrixXd computeDP2DF(int tetID);
	Eigen::Matrix3d restoreMatrix33fromTeranVector(Eigen::VectorXd v);
	const std::array<int, 9> m_matrix33fromTeran; // transfer teran's order to 3*3 matrix row major

	// compute the diagonal entries of the dP(Fhat)/dFij
	void computeDPFhat2DFij(int tetID, const double * hessian, int i, int j, double *dPFhatdFij_diagonal);

};

