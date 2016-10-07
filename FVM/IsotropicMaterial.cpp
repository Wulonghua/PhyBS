#include "IsotropicMaterial.h"


IsotropicMaterial::IsotropicMaterial() :m_eps_singularvalue(1e-8)
{
}


IsotropicMaterial::~IsotropicMaterial()
{
}

/**********see [Teran. 2004], compute F_hat and make sure U,V are real rotation matrix.********/
void IsotropicMaterial::computeSVD33modified(Eigen::Matrix3d F, Eigen::Vector3d &S, Eigen::Matrix3d &U, Eigen::Matrix3d &V)
{
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
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

void IsotropicMaterial::computeFhatsInvariants()
{
	int n = m_tetModel->getTetsNum();
	Eigen::Matrix3d F, U, V;
	Eigen::Vector3d Fhat;
	double sigma1sq, sigma2sq, sigma3sq;
	
	for (int i = 0; i < n; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i);
		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;
		m_Us.block<3, 3>(0, i * 3) = U;
		m_Vs.block<3, 3>(0, i * 3) = V;

		sigma1sq = Fhat(0)*Fhat(0);
		sigma2sq = Fhat(1)*Fhat(1);
		sigma3sq = Fhat(2)*Fhat(2);

		m_Invariants(0, i) = sigma1sq + sigma2sq + sigma3sq;
		m_Invariants(1, i) = sigma1sq * sigma1sq + sigma2sq * sigma2sq + sigma3sq * sigma3sq;
		m_Invariants(2, i) = sigma1sq * sigma2sq * sigma3sq;
	}
}