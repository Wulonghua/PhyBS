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

Eigen::MatrixXd IsotropicNeohookeanMaterial::computeDP2DF(int tetID)
{
	// first compute dP/dF at Fhat. See[Teran 05] Section 8
	Eigen::MatrixXd dPdFatFhat = Eigen::MatrixXd::Zero(9, 9);
	double invariantIII = m_Invariants(2, tetID);
	Eigen::Vector3d gradient = computeEnergy2InvariantsGradient(tetID, m_Invariants.col(tetID));
	double hessianIIIsq = computeEnergy2InvariantsHessian(tetID, m_Invariants.col(tetID))(2, 2);
	double sigma11 = m_Fhats(0, tetID) * m_Fhats(0, tetID);
	double sigma12 = m_Fhats(0, tetID) * m_Fhats(1, tetID);
	double sigma13 = m_Fhats(0, tetID) * m_Fhats(2, tetID);
	double sigma22 = m_Fhats(1, tetID) * m_Fhats(1, tetID);
	double sigma23 = m_Fhats(1, tetID) * m_Fhats(2, tetID);
	double sigma33 = m_Fhats(2, tetID) * m_Fhats(2, tetID);
	double alpha = 2.0 * gradient(0);
	double beta = -2.0 * invariantIII * gradient(2);
	double gamma = 4.0 * invariantIII * (invariantIII*hessianIIIsq + gradient(2));
	
	Eigen::Matrix3d A;
	Eigen::Matrix2d B12, B13, B23;
	A(0, 0) = alpha + (beta + gamma) / sigma11;
	A(0, 1) = A(1, 0) = gamma / sigma12;
	A(0, 2) = A(2, 0) = gamma / sigma13;
	A(1, 1) = alpha + (beta + gamma) / sigma22;
	A(1, 2) = A(2, 1) = gamma / sigma23;
	A(2, 2) = alpha + (beta + gamma) / sigma33;

	B12(0, 0) = B12(1, 1) = alpha;
	B12(0, 1) = B12(1, 0) = beta / sigma12;

	B13(0, 0) = B13(1, 1) = alpha;
	B13(0, 1) = B13(1, 0) = beta / sigma13;

	B23(0, 0) = B23(1, 1) = alpha;
	B23(0, 1) = B23(1, 0) = beta / sigma23;

	dPdFatFhat.block<3, 3>(0, 0) = A;
	dPdFatFhat.block<2, 2>(3, 3) = B12;
	dPdFatFhat.block<2, 2>(5, 5) = B13;
	dPdFatFhat.block<2, 2>(7, 7) = B23;

	// Then compute dP/dF using [Teran 05] equation (2)
	Eigen::MatrixXd dPdF = Eigen::MatrixXd::Zero(9, 9);
	Eigen::Matrix3d eij = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dPdFij = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dPdFij_middle = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dPdFij_t = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d U = m_Us.block<3, 3>(0, 3 * tetID);
	Eigen::Matrix3d V = m_Vs.block<3, 3>(0, 3 * tetID);
	Eigen::Matrix3d UT = U.transpose();
	Eigen::Matrix3d VT = V.transpose();
	Eigen::Matrix3d UTeijV;
	Eigen::Matrix3d subTensor;

	for (int fi = 0; fi < 3; ++fi)
	{
		for (int fj = 0; fj < 3; ++fj)
		{
			eij(fi, fj) = 1;
			UTeijV = UT*eij*V;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					subTensor = restoreMatrix33fromTeranVector(dPdFatFhat.row(i * 3 + j));
					dPdFij_middle += UTeijV(i,j) * subTensor;
				}
			}
			dPdFij = U * dPdFij_middle * VT;
			dPdFij_t = dPdFij.transpose();
			for (int r = 0; r < 9; ++r)
				dPdF(r, fi * 3 + fj) = dPdFij_t.data()[r];
			eij(fi, fj) = 0;
			dPdFij_middle = Eigen::Matrix3d::Zero();
		}
	}
	return dPdF;
}

Eigen::MatrixXd IsotropicNeohookeanMaterial::computeStiffnessMatrix(int tetID)
{
	Eigen::MatrixXd dPdF = computeDP2DF(tetID); // 9*9
	Eigen::Matrix3d BT = m_tetModel->getAN(tetID).transpose();
	Eigen::MatrixXd dGdF(9, 9);
	dGdF.block<3, 9>(0, 0) = BT * dPdF.block<3, 9>(0, 0);
	dGdF.block<3, 9>(3, 0) = BT * dPdF.block<3, 9>(3, 0);
	dGdF.block<3, 9>(6, 0) = BT * dPdF.block<3, 9>(6, 0);


	Eigen::Matrix3d DmInvT = m_tetModel->getDmInv(tetID).transpose();
	Eigen::Vector3d v = -1.0 * DmInvT.rowwise().sum();

	Eigen::MatrixXd dFdx = Eigen::MatrixXd::Zero(9, 9);
	dFdx.block<3, 1>(0, 0) = dFdx.block<3, 1>(3, 1) = dFdx.block<3, 1>(6, 2) = v;
	dFdx.block<3, 1>(0, 3) = dFdx.block<3, 1>(3, 4) = dFdx.block<3, 1>(6, 5) = DmInvT.col(0);
	dFdx.block<3, 1>(0, 6) = dFdx.block<3, 1>(3, 7) = dFdx.block<3, 1>(6, 8) = DmInvT.col(1);
	dFdx.block<3, 1>(0, 9) = dFdx.block<3, 1>(3, 10) = dFdx.block<3, 1>(6, 11) = DmInvT.col(2);

	Eigen::MatrixXd dGdx = dGdF *dFdx;

}

Eigen::Matrix3d IsotropicNeohookeanMaterial::restoreMatrix33fromTeranVector(Eigen::Vector3d v)
{
	Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
	double matrix33fromTeran[9] = { 0, 3, 5, 4, 1, 7, 6, 8, 2 };
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			mat(i, j) = v(i * 3 + j);
		}
	}
	return mat;
}