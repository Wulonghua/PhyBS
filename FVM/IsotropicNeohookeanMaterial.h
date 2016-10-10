#pragma once
#include "IsotropicMaterial.h"
#include <Eigen\Dense>
#include <algorithm>
#include <cmath>
#include <array>

class IsotropicNeohookeanMaterial :
	public IsotropicMaterial
{
public:
	IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh);
	virtual ~IsotropicNeohookeanMaterial();

	virtual Eigen::Vector3d computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants);
	virtual Eigen::Matrix3d computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants);

	Eigen::MatrixXd computeInnerForcesfromFhats();

	// compute stiffness Matrix
	Eigen::MatrixXd computeStiffnessMatrix(int tetID);



private:
	// compute dP/dF
	Eigen::MatrixXd computeDP2DF(int tetID);
	Eigen::Matrix3d restoreMatrix33fromTeranVector(Eigen::VectorXd v);
	const std::array<int,9> m_matrix33fromTeran; // transfer teran's order to 3*3 matrix row major
};

