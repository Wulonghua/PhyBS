#include "IsotropicNeohookeanMaterial.h"


IsotropicNeohookeanMaterial::IsotropicNeohookeanMaterial(std::shared_ptr<TetMesh> tetMesh)
	:m_matrix33fromTeran({ { 0, 3, 5, 4, 1, 7, 6, 8, 2 } }) // compromised way to initialize: VS2013 does not fully support c++11
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
	std::cout << "Isotrpic Neo-hookean Material initialized."<<std::endl;
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

Eigen::MatrixXd IsotropicNeohookeanMaterial::computeInnerForcesfromFhats()
{
	computeFhatsInvariants();
	int n = m_tetModel->getTetsNum();
	Eigen::Vector3d Phat,Fhat,Fhat_inverseT;
	double mu, lambda, I3;
	Eigen::Matrix3d P,U,V,forces;

	m_tetModel->initForcesFromGravity();
	for (int i = 0; i < n; ++i)
	{
		// compute First Piola-Kirchhoff stress based on diagonalized F.
		// Use the equation from SIGGRAPH course note[Sifaki 2012] Page 24
		Fhat = m_Fhats.col(i);
		Fhat_inverseT = Fhat.cwiseInverse();

		//std::cout << "Fhat: " << std::endl;
		//std::cout << Fhat << std::endl;
		//std::cout << "Fhat inverseT" << std::endl;
		//std::cout << Fhat_inverseT << std::endl;

		mu = m_mus[i];
		lambda = m_lambdas[i];
		I3 = m_Invariants(2, i);

		//std::cout << "I3: " << I3 << std::endl;
		Phat = (Fhat - Fhat_inverseT) * mu + 0.5 * lambda * std::log(I3) * Fhat_inverseT;

		//double tmp = Phat[0];
		//std::cout << "Phat:" << std::endl;
		//std::cout << Phat << std::endl;
		
		// equation 1 in [Teran 04] P = U * Phat * V^T
		U = m_Us.block<3, 3>(0, 3 * i);
		V = m_Vs.block<3, 3>(0, 3 * i);
		P = U * Phat.asDiagonal() * V.transpose();

		//std::cout << "U: " << std::endl;
		//std::cout << U << std::endl;

		//std::cout << "V: " << std::endl;
		//std::cout << V << std::endl;

		//std::cout << "P: " << std::endl;
		//std::cout << P << std::endl;

		forces = P * m_tetModel->getAN(i);

		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(0));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0), -(forces.rowwise().sum()));
	}

	return m_tetModel->getForces();
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

	//std::cout << dPdFatFhat << std::endl;

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
					subTensor = restoreMatrix33fromTeranVector(dPdFatFhat.row(m_matrix33fromTeran[i * 3 + j]));
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
	//m_tetModel->writeMatrix("dPdF.csv", dPdF);

	Eigen::Matrix3d BT = m_tetModel->getAN(tetID).transpose();
	Eigen::MatrixXd dGdF(9, 9);
	dGdF.block<3, 9>(0, 0) = BT * dPdF.block<3, 9>(0, 0);
	dGdF.block<3, 9>(3, 0) = BT * dPdF.block<3, 9>(3, 0);
	dGdF.block<3, 9>(6, 0) = BT * dPdF.block<3, 9>(6, 0);


	Eigen::Matrix3d DmInvT = m_tetModel->getDmInv(tetID).transpose();
	Eigen::Vector3d v = -1.0 * DmInvT.rowwise().sum();

	Eigen::MatrixXd dFdx = Eigen::MatrixXd::Zero(9, 12);
	dFdx.block<3, 1>(0, 0) = dFdx.block<3, 1>(3, 1) = dFdx.block<3, 1>(6, 2) = v;
	dFdx.block<3, 1>(0, 3) = dFdx.block<3, 1>(3, 4) = dFdx.block<3, 1>(6, 5) = DmInvT.col(0);
	dFdx.block<3, 1>(0, 6) = dFdx.block<3, 1>(3, 7) = dFdx.block<3, 1>(6, 8) = DmInvT.col(1);
	dFdx.block<3, 1>(0, 9) = dFdx.block<3, 1>(3, 10) = dFdx.block<3, 1>(6, 11) = DmInvT.col(2);

	Eigen::MatrixXd dGdx = dGdF * dFdx;

	Eigen::MatrixXd dfdx = Eigen::MatrixXd::Zero(12, 12);

	dfdx.row(0) = -dGdx.row(0) - dGdx.row(1) - dGdx.row(2);
	dfdx.row(1) = -dGdx.row(3) - dGdx.row(4) - dGdx.row(5);
	dfdx.row(2) = -dGdx.row(6) - dGdx.row(7) - dGdx.row(8);
	
	int convert_idx[9] = {0,3,6,1,4,7,2,5,8};
	for (int i = 0; i < 9; ++i)
	{
		dfdx.row(i + 3) = dGdx.row(convert_idx[i]);
	}

	// test
	//m_tetModel->writeMatrix("mat.csv", dfdx);
	return dfdx;
}

Eigen::SparseMatrix<double> IsotropicNeohookeanMaterial::computeGlobalStiffnessMatrix()
{
	int n = m_tetModel->getNodesNum();
	int m = m_tetModel->getTetsNum();
	Eigen::SparseMatrix<double> gK(3 * n, 3 * n);
	Eigen::MatrixXd K;
	int Ki, Kj, gKi,gKj;

	for (int i = 0; i < m; ++i)
	{
		K = computeStiffnessMatrix(i);
		for (int fi = 0; fi < 4; ++fi)
		{
			for (int fj = 0; fj < 3; ++fj)
			{
				Ki = fi * 3 + fj;
				gKi = m_tetModel->getNodeGlobalIDinTet(i, fi) * 3 + fj;
				for (int ni = 0; ni < 4; ++ni)
				{
					for (int nj = 0; nj < 3; ++nj)
					{
						Kj = ni * 3 + nj;
						gKj = m_tetModel->getNodeGlobalIDinTet(i, ni) * 3 + nj;
						gK.coeffRef(gKi, gKj) += K(Ki,Kj);
					}
				}
			}
		}
	}

	return gK;
}

Eigen::Matrix3d IsotropicNeohookeanMaterial::restoreMatrix33fromTeranVector(Eigen::VectorXd v)
{
	Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			mat(i, j) = v(m_matrix33fromTeran[i * 3 + j]);
		}
	}
	return mat;
}