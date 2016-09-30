#include "IsotropicMaterial.h"


IsotropicMaterial::IsotropicMaterial() :m_eps_singularvalue(1e-8)
{
}


IsotropicMaterial::~IsotropicMaterial()
{
}

void IsotropicMaterial::computeSVD33modified(Eigen::Matrix3d F, Eigen::Vector3d &S, Eigen::Matrix3d &U, Eigen::Matrix3d &V)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
	S = svd.singularValues();
	V = svd.matrixV();


	// Fix V if determinant of V is equal to -1
	if (V.determinant() < 0) // ==-1
	{
		V.col(0) = V.col(0) * (-1);
	}

	Eigen::Matrix3d ss = Eigen::Matrix3d::Zero();
	ss(0, 0) = S(0) > m_eps_singularvalue ? 1 / S(0) : 0.0;
	ss(1, 1) = S(1) > m_eps_singularvalue ? 1 / S(1) : 0.0;
	ss(2, 2) = S(2) > m_eps_singularvalue ? 1 / S(2) : 0.0;

	U = F * V * ss;
	//The returned singular values are already sorted in descending order if use Eigen lib function
	//Fix U is certain singularvalue is below epsilon or equal to 0.
	if (S(0) < m_eps_singularvalue) //all singular values are equal to 0 or below epsilon
	{
		U = Eigen::Matrix3d::Identity();
	}
	else if (S(1) < m_eps_singularvalue) // two singular values are equal to 0 or below epsilon
	{
		Eigen::Vector3d v1 = U.col(0).unitOrthogonal();
		Eigen::Vector3d v2 = U.col(0).cross(v1).normalized();
		U.col(1) = v1;
		U.col(2) = v2;
	}
	else if (S(2) < m_eps_singularvalue) // one singular value is equal to 0 or below epsilon
	{
		U.col(2) = U.col(0).cross(U.col(1)).normalized();
	}

	//Fix U and negate minimal singularvalue if dterminant of U is equal to -1
	if (U.determinant() < 0) // ==-1
	{
		U.col(2) = U.col(2) * (-1);
		S(2) *= -1;
	}
}