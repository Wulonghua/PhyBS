#include "IsotropicNeohookeanMaterial.h"


IsotropicNeohookeanMaterial::IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh)
{
	int n_tets = tetMesh->getTetsNum();
	m_mus.resize(n_tets);
	m_lambdas.resize(n_tets);
	m_elasticEnergys.resize(n_tets);
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

void IsotropicNeohookeanMaterial::computeEnergy2FhatGradient(int tetID, const float *Fhats, float *gradient)
{
	float mu = m_mus[tetID];
	float lambda = m_lambdas[tetID];
	float tmp = lambda *  (std::log(Fhats[0]) + std::log(Fhats[1]) + std::log(Fhats[2])) - mu;
	//std::cout << "mu: " << mu << std::endl;
	//std::cout << "lambda: " << lambda << std::endl;

	for (int i = 0; i < 3; ++i)
		gradient[i] = mu*Fhats[i] + tmp / Fhats[i];
}

// *hessian is the one dimentional representation of the row-major 3*3 hessian matrix
void IsotropicNeohookeanMaterial::computeEnergy2FhatHessian(int tetID, const float *Fhats, float *hessian)
{
	float mu = m_mus[tetID];
	float lambda = m_lambdas[tetID];
	
	float tmp  =  lambda * (std::log(Fhats[0]) + std::log(Fhats[1]) + std::log(Fhats[2]));
	float tmp1 = mu + lambda - tmp;
	float tmp2 = tmp - mu;

	float inv_lambda1 = 1 / Fhats[0];
	float inv_lambda2 = 1 / Fhats[1];
	float inv_lambda3 = 1 / Fhats[2];

		// hessian(1,1)
	hessian[0] = mu + tmp1 * inv_lambda1 * inv_lambda1;
		// hessian(2,2)
	hessian[4] = mu + tmp1 * inv_lambda2 * inv_lambda2;
		// hessian(3,3)
	hessian[8] = mu + tmp1 * inv_lambda3 * inv_lambda3;

		// hessian(1,2) = hessian(2,1)
	hessian[1] = hessian[3] = (tmp1 + tmp2) * inv_lambda1 * inv_lambda2;
	
		// hessian(1,3) = hessian(3,1)
	hessian[2] = hessian[6] = (tmp1 + tmp2) * inv_lambda1 * inv_lambda3;
	
		// hessian(2,3) = hessian(3,2)
	hessian[5] = hessian[7] = (tmp1 + tmp2) * inv_lambda2 * inv_lambda3;
}

Eigen::MatrixXf IsotropicNeohookeanMaterial::computeInnerForceFromPos(const Eigen::MatrixXf & pos)
{
	Eigen::Matrix3f F,F_invT,P,forces;
	int n_tets = m_tetModel->getTetsNum();

	m_tetModel->initForcesFromGravityExternals();
	for (int i = 0; i < n_tets; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i, pos);
		F_invT = F.inverse().transpose();

		P = (F - F_invT)*m_mus[i] + m_lambdas[i] * std::log(F.determinant())*F_invT;
		forces = P * m_tetModel->getAN(i);

		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(0));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0), -(forces.rowwise().sum()));
	}
	m_tetModel->addFixNodeSpringForces();
	return m_tetModel->getForces();
}

void IsotropicNeohookeanMaterial::computeElasticEnergyFromPos(const Eigen::MatrixXf & pos)
{
	Eigen::Matrix3f F;
	float I1, logJ;
	int n_tets = m_tetModel->getTetsNum();
	for (int i = 0; i < n_tets; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i, pos);
		I1 = F.squaredNorm();
		logJ = std::log(F.determinant());
		m_elasticEnergys[i] = 0.5 * m_mus[i] * (I1 - 3) - m_mus[i] * logJ + 0.5 * m_lambdas[i] * logJ * logJ;
	}
}


