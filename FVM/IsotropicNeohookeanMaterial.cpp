#include "IsotropicNeohookeanMaterial.h"


IsotropicNeohookeanMaterial::IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh)
{
	int n_tets = tetMesh->getTetsNum();
	m_mus.resize(n_tets);
	m_lambdas.resize(n_tets);
	m_Fhats = Eigen::MatrixXf::Zero(3, n_tets);
	m_Us = Eigen::MatrixXf::Zero(3, 3 * n_tets);
	m_Vs = Eigen::MatrixXf::Zero(3, 3 * n_tets);
	m_Invariants = Eigen::MatrixXf::Zero(3, n_tets);

	float mu = 0.5 * tetMesh->getE() / (1 + tetMesh->getNu());
	float lambda = mu * 2 / (1 - 2 * tetMesh->getNu());

	//currently only homogeneous material is used. 
	std::fill(m_mus.begin(),m_mus.end(),mu);
	std::fill(m_lambdas.begin(), m_lambdas.end(), lambda);

	m_tetModel = tetMesh;
	allocateGlobalStiffnessMatrix();
	std::cout << "Isotrpic Neo-hookean Material initialized."<<std::endl;
}


IsotropicNeohookeanMaterial::~IsotropicNeohookeanMaterial()
{
}

Eigen::Vector3f IsotropicNeohookeanMaterial::computeEnergy2InvariantsGradient(int tetID, Eigen::Vector3f invariants)
{
	// I_1: invariants(0);	I_2: invariants(1);	I_3: invariants(2)
	Eigen::Vector3f gradient = Eigen::Vector3f::Zero();
	gradient(0) = m_mus[tetID] * 0.5;
	gradient(2) = -m_mus[tetID] * 0.5 / invariants(2) + 0.25 * m_lambdas[tetID] * std::log(invariants(2)) / invariants(2);
	return gradient;
}

Eigen::Matrix3f IsotropicNeohookeanMaterial::computeEnergy2InvariantsHessian(int tetID, Eigen::Vector3f invariants)
{
	Eigen::Matrix3f hessian = Eigen::Matrix3f::Zero();
	float invI3sq = 1 / (invariants(2) * invariants(2));
	hessian(2, 2) = 0.5 * m_mus[tetID] * invI3sq + 0.25 * m_lambdas[tetID] * invI3sq * (1 - std::log(invariants(2)));

	return hessian;
}



