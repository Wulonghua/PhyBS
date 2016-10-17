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

protected:
	virtual double fEnergy(double x, const double & mu, const double & lambda) { return 0.5*mu *(x*x - 1); }
	virtual double gEnergy(double x, const double & mu, const double & lambda) { return 0; }
	virtual double hEnergy(double x, const double & mu, const double & lambda)
	{
		double lgx = std::log(x);
		return 0.5 * lambda * (lgx*lgx) - mu*lgx;
	}
	virtual double dfEnergy(double x, const double & mu, const double & lambda) { return mu * x; }
	virtual double dgEnergy(double x, const double & mu, const double & lambda) { return 0; }
	virtual double dhEnergy(double x, const double & mu, const double & lambda) { return (lambda*std::log(x) - mu) / x; }
	virtual double ddfEnergy(double x, const double & mu, const double & lambda) { return mu; }
	virtual double ddgEnergy(double x, const double & mu, const double & lambda) { return 0; }
	virtual double ddhEnergy(double x, const double & mu, const double & lambda) { return (mu + lambda - lambda*std::log(x)) / (x*x); }
};

