#pragma once
#include <Eigen\Dense>

class IsotropicMaterial
{
public:
	IsotropicMaterial();
	virtual ~IsotropicMaterial();

	//virtual double computeEnergy(int tetID, Eigen::Vector3d invariants);
	virtual Eigen::Vector3d computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants)=0;
	virtual Eigen::Matrix3d computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants)=0;
};

