#pragma once
#include "IsotropicMaterial.h"
#include <Eigen\Dense>
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
	Eigen::MatrixXd computeDP2DF(int tetID);
	void computeInnerForcesfromFhats();

private:
	Eigen::Matrix3d restoreMatrix33fromTeranVector(Eigen::Vector3d v);

};

