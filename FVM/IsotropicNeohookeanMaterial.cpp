#include "IsotropicNeohookeanMaterial.h"


IsotropicNeohookeanMaterial::IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh)
{
	int n_tets = tetMesh->getTetsNum();
	m_mus.resize(n_tets);
	m_lambdas.resize(n_tets);
	m_Fhats = Eigen::MatrixXd::Zero(3, n_tets);
	m_Us = Eigen::MatrixXd::Zero(3, 3 * n_tets);
	m_Vs = Eigen::MatrixXd::Zero(3, 3 * n_tets);
	m_Invariants = Eigen::MatrixXd::Zero(3, n_tets);

	double mu = 0.5 * tetMesh->getE() / (1 + tetMesh->getNu());
	double lambda = mu * 2 / (1 - 2 * tetMesh->getNu());

	//currently only homogeneous material is used. 
	std::fill(m_mus.begin(),m_mus.end(),mu);
	std::fill(m_lambdas.begin(), m_lambdas.end(), lambda);

	m_tetModel = tetMesh;
}


IsotropicNeohookeanMaterial::~IsotropicNeohookeanMaterial()
{
}

Eigen::Vector3d IsotropicNeohookeanMaterial::computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3d invariants)
{
	// I_1: invariants(0);	I_2: invariants(1);	I_3: invariants(2)
	Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
	gradient(0) = m_mus[tetID] * 0.5;
	gradient(2) = -m_mus[tetID] * 0.5 / invariants(2) + 0.25 * m_lambdas[tetID] * std::log(invariants(2)) / invariants(2);
	return gradient;
}

Eigen::Matrix3d IsotropicNeohookeanMaterial::computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3d invariants)
{
	Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();
	double invI3sq = 1 / (invariants(2) * invariants(2));
	hessian(2, 2) = 0.5 * m_mus[tetID] * invI3sq + 0.25 * m_lambdas[tetID] * invI3sq * (1 - std::log(invariants(2)));

	return hessian;
}

void IsotropicNeohookeanMaterial::computeInnerForcesfromFhats()
{
	computeFhatsInvariants();
	int n = m_tetModel->getTetsNum();
	Eigen::Vector3d Phat,Fhat,F_inverseT;
	double mu, lambda, J;
	Eigen::Matrix3d P,U,V,forces;

	for (int i = 0; i < n; ++i)
	{
		// compute First Piola-Kirchhoff stress based on diagonalized F.
		// Use the equation from SIGGRAPH course note[Sifaki 2012] Page 24
		Fhat = m_Fhats.col(i);
		F_inverseT = Fhat.cwiseInverse();
		mu = m_mus[i];
		lambda = m_lambdas[i];
		J = std::sqrt(m_Invariants(2, i));
		Phat = mu * (Fhat - mu * F_inverseT) + lambda * std::log(J) * F_inverseT;
		
		// equation 1 in [Teran 04] P = U * Fhat * V^T
		U = m_Us.block<3, 3>(0, 3 * i);
		V = m_Vs.block<3, 3>(0, 3 * i);
		P.col(0) = U.col(0) * Phat(0);
		P.col(1) = U.col(1) * Phat(1);
		P.col(2) = U.col(2) * Phat(2);
		P = P * V.transpose();

		forces = P * m_tetModel->getAN(i);

		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(3));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0),
								 -(forces.col(0) + forces.col(1) + forces.col(2)));	
	}
}

Eigen::MatrixXd IsotropicNeohookeanMaterial::computePDP2PDF(int tetID)
{

	return Eigen::Matrix3d::Zero();
}