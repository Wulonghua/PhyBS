#pragma once
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen\Dense>
#include <Eigen\Core>
#include <Eigen\svd>
#include <vector>
#include <memory>

#include "TetMesh.h"

class IsotropicMaterial
{
public:
	IsotropicMaterial();
	virtual ~IsotropicMaterial();

	//virtual double computeEnergy(int tetID, Eigen::Vector3d invariants);
	virtual Eigen::Vector3d computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants)=0;
	virtual Eigen::Matrix3d computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants)=0;

	//see [Teran. 2004], make sure U,V are real rotation matrix.
	void computeSVD33modified(Eigen::Matrix3d F, Eigen::Vector3d &S, Eigen::Matrix3d &U, Eigen::Matrix3d &V);

	const double m_eps_singularvalue;

protected:
	std::vector<double> m_mus;
	std::vector<double> m_lambdas;
	std::shared_ptr<TetMesh> m_tetModel;
};

