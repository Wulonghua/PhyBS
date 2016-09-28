#pragma once
#include "IsotropicMaterial.h"
#include "TetMesh.h"
#include <Eigen\Dense>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>

class IsotropicNeohookeanMaterial :
	public IsotropicMaterial
{
public:
	IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh);
	virtual ~IsotropicNeohookeanMaterial();

	virtual Eigen::Vector3d computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants);
	virtual Eigen::Matrix3d computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants);
	// compute dP/dF
	Eigen::MatrixXd computePDP2PDF(int tetID);

private:
	std::vector<double> m_mus;
	std::vector<double> m_lambdas;
	std::shared_ptr<TetMesh> m_tetModel;

};

