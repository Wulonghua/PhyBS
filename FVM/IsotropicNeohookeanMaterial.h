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

	virtual Eigen::Vector3f computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3f invariants);
	virtual Eigen::Matrix3f computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3f invariants);

protected:
	virtual float fEnergy(float x, const float & mu, const float & lambda) { return 0.5*mu *(x*x - 1); }
	virtual float gEnergy(float x, const float & mu, const float & lambda) { return 0; }
	virtual float hEnergy(float x, const float & mu, const float & lambda)
	{
		float lgx = std::log(x);
		return 0.5 * lambda * (lgx*lgx) - mu*lgx;
	}
	virtual float dfEnergy(float x, const float & mu, const float & lambda) { return mu * x; }
	virtual float dgEnergy(float x, const float & mu, const float & lambda) { return 0; }
	virtual float dhEnergy(float x, const float & mu, const float & lambda) { return (lambda*std::log(x) - mu) / x; }
	virtual float ddfEnergy(float x, const float & mu, const float & lambda) { return mu; }
	virtual float ddgEnergy(float x, const float & mu, const float & lambda) { return 0; }
	virtual float ddhEnergy(float x, const float & mu, const float & lambda) { return (mu + lambda - lambda*std::log(x)) / (x*x); }
};

