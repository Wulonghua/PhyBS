#include "IsotropicMaterial.h"

IsotropicMaterial::IsotropicMaterial() :
m_eps_singularvalue(1e-8), 
m_matrix33fromTeran({ { 0, 3, 5, 4, 1, 7, 6, 8, 2 } }) // compromised way to initialize: VS2013 does not fully support c++11
{
	//m_timeTest.start();
}


IsotropicMaterial::~IsotropicMaterial()
{
}

/**********see [Teran. 2004], compute F_hat and make sure U,V are real rotation matrix.********/
void IsotropicMaterial::computeSVD33modified(Eigen::Matrix3f F, Eigen::Vector3f &S, Eigen::Matrix3f &U, Eigen::Matrix3f &V)
{
	Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	S = svd.singularValues();
	//for (int i = 0; i < 3; ++i)
	//	S(i) = m_tetModel->fixPrecision(S(i));
	V = svd.matrixV();

	// Fix V if determinant of V is equal to -1
	if (V.determinant() < 0) // ==-1
	{
		V.col(0) = V.col(0) * (-1);
	}

	Eigen::Matrix3f ss = Eigen::Matrix3f::Zero();
	ss(0, 0) = S(0) > m_eps_singularvalue ? 1 / S(0) : 0.0;
	ss(1, 1) = S(1) > m_eps_singularvalue ? 1 / S(1) : 0.0;
	ss(2, 2) = S(2) > m_eps_singularvalue ? 1 / S(2) : 0.0;

	U = F * V * ss;
	//The returned singular values are already sorted in descending order if use Eigen lib function
	//Fix U if certain singularvalue is below epsilon or equal to 0.
	if (S(0) < m_eps_singularvalue) //all singular values are equal to 0 or below epsilon
	{
		U = Eigen::Matrix3f::Identity();
	}
	else if (S(1) < m_eps_singularvalue) // two singular values are equal to 0 or below epsilon
	{
		Eigen::Vector3f v1 = U.col(0).unitOrthogonal();
		Eigen::Vector3f v2 = U.col(0).cross(v1).normalized();
		U.col(1) = v1;
		U.col(2) = v2;
	}
	else if (S(2) < m_eps_singularvalue) // one singular value is equal to 0 or below epsilon
	{
		U.col(2) = U.col(0).cross(U.col(1)).normalized();
	}

	//Fix U and negate minimal singularvalue if determinant of U is equal to -1
	if (U.determinant() < 0) // ==-1
	{
		U.col(2) = U.col(2) * (-1);
		S(2) *= -1;
	}
}

void IsotropicMaterial::computeFhats()
{
	int n = m_tetModel->getTetsNum();
	Eigen::Matrix3f F, U, V;
	Eigen::Vector3f Fhat;

	for (int i = 0; i < n; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i);

		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;
		//float t = Fhat[0];

		m_Us.block<3, 3>(0, i * 3) = U;
		m_Vs.block<3, 3>(0, i * 3) = V;
	}
}

void IsotropicMaterial::computeFhats(int num_Threads)
{
	int n = m_tetModel->getTetsNum();

	omp_set_num_threads(num_Threads);
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		Eigen::Matrix3f F, U, V;
		Eigen::Vector3f Fhat;

		F = m_tetModel->computeDeformationGradient(i);
		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;

		m_Us.block<3, 3>(0, i * 3) = U;
		m_Vs.block<3, 3>(0, i * 3) = V;
	}
}


void IsotropicMaterial::computeFhatsInvariants()
{
	int n = m_tetModel->getTetsNum();
	Eigen::Matrix3f F, U, V;
	Eigen::Vector3f Fhat;
	float sigma1sq, sigma2sq, sigma3sq;
	
	for (int i = 0; i < n; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i);

		//std::cout << "F: " << std::endl;
		//std::cout << F << std::endl;

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

void IsotropicMaterial::computeFhatsInvariants(int num_Threads)
{
	int n = m_tetModel->getTetsNum();

	omp_set_num_threads(num_Threads);
	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		Eigen::Matrix3f F, U, V;
		Eigen::Vector3f Fhat;

		F = m_tetModel->computeDeformationGradient(i);

		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;
		float t = Fhat[0];

		m_Us.block<3, 3>(0, i * 3) = U;
		m_Vs.block<3, 3>(0, i * 3) = V;

		float sigma1sq = Fhat(0)*Fhat(0);
		float sigma2sq = Fhat(1)*Fhat(1);
		float sigma3sq = Fhat(2)*Fhat(2);

		m_Invariants(0, i) = sigma1sq + sigma2sq + sigma3sq;
		m_Invariants(1, i) = sigma1sq * sigma1sq + sigma2sq * sigma2sq + sigma3sq * sigma3sq;
		m_Invariants(2, i) = sigma1sq * sigma2sq * sigma3sq;
	}
}

Eigen::MatrixXf IsotropicMaterial::computeDP2DF(int tetID)
{
	// first compute dP/dF at Fhat. See[Teran 05] Section 8
	Eigen::MatrixXf dPdFatFhat = Eigen::MatrixXf::Zero(9, 9);
	float invariantIII = m_Invariants(2, tetID);
	Eigen::Vector3f gradient = computeEnergy2InvariantsGradient(tetID, m_Invariants.col(tetID));
	float hessianIIIsq = computeEnergy2InvariantsHessian(tetID, m_Invariants.col(tetID))(2, 2);
	float sigma11 = m_Fhats(0, tetID) * m_Fhats(0, tetID);
	float sigma12 = m_Fhats(0, tetID) * m_Fhats(1, tetID);
	float sigma13 = m_Fhats(0, tetID) * m_Fhats(2, tetID);
	float sigma22 = m_Fhats(1, tetID) * m_Fhats(1, tetID);
	float sigma23 = m_Fhats(1, tetID) * m_Fhats(2, tetID);
	float sigma33 = m_Fhats(2, tetID) * m_Fhats(2, tetID);
	float alpha = 2.0 * gradient(0);
	float beta = -2.0 * invariantIII * gradient(2);
	float gamma = 4.0 * invariantIII * (invariantIII*hessianIIIsq + gradient(2));

	Eigen::Matrix3f A;
	Eigen::Matrix2f B12, B13, B23;
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
	Eigen::MatrixXf dPdF = Eigen::MatrixXf::Zero(9, 9);
	Eigen::Matrix3f eij = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f dPdFij = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f dPdFij_middle = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f dPdFij_t = Eigen::Matrix3f::Zero();
	Eigen::Matrix3f U = m_Us.block<3, 3>(0, 3 * tetID);
	Eigen::Matrix3f V = m_Vs.block<3, 3>(0, 3 * tetID);
	Eigen::Matrix3f UT = U.transpose();
	Eigen::Matrix3f VT = V.transpose();
	Eigen::Matrix3f UTeijV;
	Eigen::Matrix3f subTensor;

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
					dPdFij_middle += UTeijV(i, j) * subTensor;
				}
			}
			dPdFij = U * dPdFij_middle * VT;
			dPdFij_t = dPdFij.transpose();
			for (int r = 0; r < 9; ++r)
				dPdF(r, fi * 3 + fj) = dPdFij_t.data()[r];
			eij(fi, fj) = 0;
			dPdFij_middle = Eigen::Matrix3f::Zero();
		}
	}
	return dPdF;
}

Eigen::MatrixXf IsotropicMaterial::computeInnerForcesfromFhats()
{
	computeFhatsInvariants();
	//computeFhats();
	int n = m_tetModel->getTetsNum();
	Eigen::Vector3f Phat, Fhat, Fhat_inverseT;
	float mu, lambda, I3;
	Eigen::Matrix3f P, U, V, forces;

	//float Phat_[3];

	m_tetModel->initForcesFromGravityExternals();
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
		//computeEnergy2FhatGradient(i, Fhat.data(), Phat_);

		//if (i == 0)
		//{
		//	std::cout << "I3: " << I3 << std::endl;
		//	std::cout << "Fhats: " << Fhat[0] << " " << Fhat[1] << " " << Fhat[2] << std::endl;
		//	std::cout << "Fhats_inverse: " << Fhat_inverseT[0] << " " << Fhat_inverseT[1] << " " << Fhat_inverseT[2] << std::endl;
		//	std::cout << "Phat: " << std::endl;

		//	for (int j = 0; j < 3; ++j)
		//		std::cout << Phat[j] << " " << Phat_[j] << std::endl;
		//}


		//float tmp = Phat[0];
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

Eigen::MatrixXf IsotropicMaterial::computeInnerForcesfromFhats(int num_Threads)
{
	computeFhatsInvariants(num_Threads);
	int n = m_tetModel->getTetsNum();
	m_tetModel->initForcesFromGravityExternals();
	
	omp_lock_t lck;
	omp_init_lock(&lck);

	omp_set_num_threads(num_Threads);
	#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		Eigen::Vector3f Phat, Fhat, Fhat_inverseT;
		float mu, lambda, I3;
		Eigen::Matrix3f P, U, V, forces;
		// compute First Piola-Kirchhoff stress based on diagonalized F.
		// Use the equation from SIGGRAPH course note[Sifaki 2012] Page 24
		Fhat = m_Fhats.col(i);
		Fhat_inverseT = Fhat.cwiseInverse();

		mu = m_mus[i];
		lambda = m_lambdas[i];
		I3 = m_Invariants(2, i);

		Phat = (Fhat - Fhat_inverseT) * mu + 0.5 * lambda * std::log(I3) * Fhat_inverseT;

		// equation 1 in [Teran 04] P = U * Phat * V^T
		U = m_Us.block<3, 3>(0, 3 * i);
		V = m_Vs.block<3, 3>(0, 3 * i);
		P = U * Phat.asDiagonal() * V.transpose();

		forces = P * m_tetModel->getAN(i);

		omp_set_lock(&lck);
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(0));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0), -(forces.rowwise().sum()));
		omp_unset_lock(&lck);
	}
	return m_tetModel->getForces();
}

Eigen::MatrixXf IsotropicMaterial::computeInnerForcesfromFhats2()
{
	computeFhats();
	int n = m_tetModel->getTetsNum();
	Eigen::Vector3f Phat, Fhat;
	float Phat_[3];

	Eigen::Matrix3f P, U, V, forces;

	m_tetModel->initForcesFromGravityExternals();
	for (int i = 0; i < n; ++i)
	{
		// compute First Piola-Kirchhoff stress based on diagonalized F.
		// Use the equation from paper [Xu et al. 2015] equation 12
		Fhat = m_Fhats.col(i);
		//Fhat_inverseT = Fhat.cwiseInverse();
		computeEnergy2FhatGradient(i, Fhat.data(), Phat_);
		Phat[0] = Phat_[0];
		Phat[1] = Phat_[1];
		Phat[2] = Phat_[2];

		// equation 1 in [Teran 04] P = U * Phat * V^T
		U = m_Us.block<3, 3>(0, 3 * i);
		V = m_Vs.block<3, 3>(0, 3 * i);
		P = U * Phat.asDiagonal() * V.transpose();

		forces = P * m_tetModel->getAN(i);

		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(0));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0), -(forces.rowwise().sum()));
	}

	return m_tetModel->getForces();
}

Eigen::MatrixXf IsotropicMaterial::computeInnerForcesfromFhats2(int num_Threads)
{
	computeFhats(num_Threads);
	int n = m_tetModel->getTetsNum();
	m_tetModel->initForcesFromGravityExternals();
	
	omp_lock_t lck;
	omp_init_lock(&lck);
	omp_set_num_threads(num_Threads);
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		Eigen::Vector3f Phat, Fhat, Fhat_inverseT;
		float Phat_[3];
		Eigen::Matrix3f P, U, V, forces;
		// compute First Piola-Kirchhoff stress based on diagonalized F.
		// Use the equation from paper [Xu et al. 2015] equation 12
		Fhat = m_Fhats.col(i);
		//Fhat_inverseT = Fhat.cwiseInverse();
		computeEnergy2FhatGradient(i, Fhat.data(), Phat_);
		Phat[0] = Phat_[0];
		Phat[1] = Phat_[1];
		Phat[2] = Phat_[2];

		// equation 1 in [Teran 04] P = U * Phat * V^T
		U = m_Us.block<3, 3>(0, 3 * i);
		V = m_Vs.block<3, 3>(0, 3 * i);
		P = U * Phat.asDiagonal() * V.transpose();

		forces = P * m_tetModel->getAN(i);

		omp_set_lock(&lck);
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 1), forces.col(0));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 2), forces.col(1));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 3), forces.col(2));
		m_tetModel->addNodeForce(m_tetModel->getNodeGlobalIDinTet(i, 0), -(forces.rowwise().sum()));
		omp_unset_lock(&lck);
	}

	return m_tetModel->getForces();
}

Eigen::MatrixXf IsotropicMaterial::computeStiffnessMatrix(int tetID)
{

	// Teran's method
	//Eigen::MatrixXf dPdF = computeDP2DF(tetID); // 9*9
	
	// Barbic's method
	float dPdF_buf[81] = { 0.0 }; //a raw array of float
	Eigen::Map<Eigen::Matrix<float, 9, 9>> dPdF(dPdF_buf);
	Eigen::Matrix3f U = m_Us.block<3, 3>(0, tetID * 3);
	Eigen::Matrix3f V = m_Vs.block<3, 3>(0, tetID * 3);
	Eigen::Vector3f Fhats = m_Fhats.col(tetID);
	computeDP2DF(tetID, U.transpose().data(), Fhats.data(), V.transpose().data(), dPdF_buf);

	Eigen::Matrix3f BT = m_tetModel->getAN(tetID).transpose();
	Eigen::MatrixXf dGdF(9, 9);
	dGdF.block<3, 9>(0, 0) = BT * dPdF.block<3, 9>(0, 0);
	dGdF.block<3, 9>(3, 0) = BT * dPdF.block<3, 9>(3, 0);
	dGdF.block<3, 9>(6, 0) = BT * dPdF.block<3, 9>(6, 0);

	Eigen::Matrix3f DmInvT = m_tetModel->getDmInv(tetID).transpose();
	Eigen::Vector3f v = -1.0 * DmInvT.rowwise().sum();

	Eigen::MatrixXf dFdx = Eigen::MatrixXf::Zero(9, 12);
	dFdx.block<3, 1>(0, 0) = dFdx.block<3, 1>(3, 1) = dFdx.block<3, 1>(6, 2) = v;
	dFdx.block<3, 1>(0, 3) = dFdx.block<3, 1>(3, 4) = dFdx.block<3, 1>(6, 5) = DmInvT.col(0);
	dFdx.block<3, 1>(0, 6) = dFdx.block<3, 1>(3, 7) = dFdx.block<3, 1>(6, 8) = DmInvT.col(1);
	dFdx.block<3, 1>(0, 9) = dFdx.block<3, 1>(3, 10) = dFdx.block<3, 1>(6, 11) = DmInvT.col(2);

	Eigen::MatrixXf dGdx = dGdF * dFdx;

	Eigen::MatrixXf dfdx = Eigen::MatrixXf::Zero(12, 12);

	dfdx.row(0) = -dGdx.row(0) - dGdx.row(1) - dGdx.row(2);
	dfdx.row(1) = -dGdx.row(3) - dGdx.row(4) - dGdx.row(5);
	dfdx.row(2) = -dGdx.row(6) - dGdx.row(7) - dGdx.row(8);

	int convert_idx[9] = { 0, 3, 6, 1, 4, 7, 2, 5, 8 };
	for (int i = 0; i < 9; ++i)
	{
		dfdx.row(i + 3) = dGdx.row(convert_idx[i]);
	}

	// test
	//m_tetModel->writeMatrix("mat.csv", dfdx);
	return -dfdx;
}

void IsotropicMaterial::allocateGlobalStiffnessMatrix()
{
	int n = m_tetModel->getNodesNum();
	int m = m_tetModel->getTetsNum();
	int gKi, gKj;

	m_globalK.resize(3*n,3*n);
	//if (n >= 10)
	//	gK.reserve(Eigen::VectorXf::Constant(3 * n, 120));
	m_reserveSize.reserve(3 * n);

	for (int i = 0; i < m; ++i)
	{

		for (int fi = 0; fi < 4; ++fi)
		{
			for (int fj = 0; fj < 3; ++fj)
			{
				gKi = m_tetModel->getNodeGlobalIDinTet(i, fi) * 3 + fj;
				for (int ni = 0; ni < 4; ++ni)
				{
					for (int nj = 0; nj < 3; ++nj)
					{
						gKj = m_tetModel->getNodeGlobalIDinTet(i, ni) * 3 + nj;
						m_globalK.coeffRef(gKi, gKj) += 1;
					}
				}
			}
		}
	}

	for (int i = 0; i < 3 * n; ++i)
		m_reserveSize.push_back(m_globalK.col(i).sum());
}

Eigen::SparseMatrix<float, Eigen::RowMajor> IsotropicMaterial::computeGlobalStiffnessMatrix()
{
	int n = m_tetModel->getNodesNum();
	int m = m_tetModel->getTetsNum();

	//Eigen::SparseMatrix<float, Eigen::RowMajor> gK(3 * n, 3 * n);
	//gK.reserve(m_reserveSize);
	Eigen::SparseMatrix<float, Eigen::RowMajor> gK = m_globalK;
	Eigen::MatrixXf K;
	int Ki, Kj, gKi, gKj;
	//m_timeTest.restart();

	//float t_K=0;
	//float t_globalK=0;

	for (int i = 0; i < m; ++i)
	{
		//m_timeTest.restart();
		K = computeStiffnessMatrix(i);
		//t_K += m_timeTest.restart();
		
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
						gK.coeffRef(gKi, gKj) += K(Ki, Kj)-1;
					}
				}
			}
		}
		//t_globalK += m_timeTest.restart();
	}

	//std::cout << "time to compute all elements' K: " << t_K << std::endl;
	//std::cout << "time to construct global K: " << t_globalK << std::endl;
	return gK;
}

//OpenMP parallel version
Eigen::SparseMatrix<float, Eigen::RowMajor> IsotropicMaterial::computeGlobalStiffnessMatrix(int num_Threads)
{
	int n = m_tetModel->getNodesNum();
	int m = m_tetModel->getTetsNum();
	Eigen::SparseMatrix<float, Eigen::RowMajor> gK = m_globalK;

	//omp_lock_t lck;
	//omp_init_lock(&lck);

	omp_set_num_threads(num_Threads);
#pragma omp parallel for
	for (int i = 0; i < m; ++i)
	{
		Eigen::MatrixXf K;
		int Ki, Kj, gKi, gKj;

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
						//omp_set_lock(&lck);
						#pragma omp atomic
						// trick here: why minus 1? because to reserve the space for non-zero entry, 
						// I added 1 to each entry when creating the initial global stiffness matrix.
						gK.coeffRef(gKi, gKj) += K(Ki, Kj) - 1;
						//omp_unset_lock(&lck);
					}
				}
			}
		}
	}
	return gK;
}

Eigen::Matrix3f IsotropicMaterial::restoreMatrix33fromTeranVector(Eigen::VectorXf v)
{
	Eigen::Matrix3f mat = Eigen::Matrix3f::Zero();
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			mat(i, j) = v(m_matrix33fromTeran[i * 3 + j]);
		}
	}
	return mat;
}

// To speed up and handle low precision of float type, 
// I manually simplify the equation rather than call d*Energy and dd*Engergy multiple times.
//void IsotropicMaterial::computeEnergy2FhatGradient(int tetID, const float *Fhats, float *gradient)
//{
//	float lame_mu = m_mus[tetID];
//	float lame_lambda = m_lambdas[tetID];
//	float lambda1 = Fhats[0];
//	float lambda2 = Fhats[1];
//	float lambda3 = Fhats[2];
//
//	float dg12 = dgEnergy(lambda1*lambda2, lame_mu, lame_lambda);
//	float dg23 = dgEnergy(lambda2*lambda3, lame_mu, lame_lambda);
//	float dg31 = dgEnergy(lambda3*lambda1, lame_mu, lame_lambda);
//
//	float dh123 = dhEnergy(lambda1*lambda2*lambda3, lame_mu, lame_lambda);
//
//	gradient[0] = dfEnergy(lambda1, lame_mu, lame_lambda) + dg12 * lambda2 + dg31 * lambda3 + dh123*lambda2*lambda3;
//	gradient[1] = dfEnergy(lambda2, lame_mu, lame_lambda) + dg23 * lambda3 + dg12 * lambda1 + dh123*lambda3*lambda1;
//	gradient[2] = dfEnergy(lambda3, lame_mu, lame_lambda) + dg31 * lambda1 + dg23 * lambda2 + dh123*lambda1*lambda2;
//}


// *hessian is the one dimentional representation of the row-major 3*3 hessian matrix
//void IsotropicMaterial::computeEnergy2FhatHessian(int tetID, const float *Fhats, float *hessian)
//{
//	float lame_mu = m_mus[tetID];
//	float lame_lambda = m_lambdas[tetID];
//
//	float lambda12 = Fhats[0] * Fhats[1];
//	float lambda23 = Fhats[1] * Fhats[2];
//	float lambda31 = Fhats[2] * Fhats[0];
//	float lambda123 = lambda12 * Fhats[2];
//
//	float dg12 = dgEnergy(lambda12, lame_mu, lame_lambda);
//	float dg23 = dgEnergy(lambda23, lame_mu, lame_lambda);
//	float dg31 = dgEnergy(lambda31, lame_mu, lame_lambda);
//
//	float ddg12 = ddgEnergy(lambda12, lame_mu, lame_lambda);
//	float ddg23 = ddgEnergy(lambda23, lame_mu, lame_lambda);
//	float ddg31 = ddgEnergy(lambda31, lame_mu, lame_lambda);
//	float dh123 = dhEnergy(lambda123, lame_mu, lame_lambda);
//	float ddh123 = ddhEnergy(lambda123, lame_mu, lame_lambda);
//
//	// hessian(1,1)
//	hessian[0] = ddfEnergy(Fhats[0], lame_mu, lame_lambda) + ddg12 * Fhats[1] * Fhats[1]
//														   + ddg31 * Fhats[2] * Fhats[2]
//														   + ddh123 * lambda23 * lambda23;
//	// hessian(2,2)
//	hessian[4] = ddfEnergy(Fhats[1], lame_mu, lame_lambda) + ddg23 * Fhats[2] * Fhats[2]
//														   + ddg12 * Fhats[0] * Fhats[0]
//														   + ddh123 * lambda31 * lambda31;
//	// hessian(3,3)
//	hessian[8] = ddfEnergy(Fhats[2], lame_mu, lame_lambda) + ddg31 * Fhats[0] * Fhats[0]
//														   + ddg23 * Fhats[1] * Fhats[1]
//														   + ddh123 * lambda12 * lambda12;
//	// hessian(1,2) = hessian(2,1)
//	hessian[1] = hessian[3] = ddg12 * lambda12 + dg12 + ddh123 * lambda23 * lambda31 + dh123 * Fhats[2];
//
//	// hessian(1,3) = hessian(3,1)
//	hessian[2] = hessian[6] = ddg31 * lambda31 + dg31 + ddh123 * lambda12 * lambda23 + dh123 * Fhats[1];
//
//	// hessian(2,3) = hessian(3,2)
//	hessian[5] = hessian[7] = ddg23 * lambda23 + dg23 + ddh123 * lambda12 * lambda31 + dh123 * Fhats[0];
//}

// see [Xu et al. 2015] Section 3.1 euqation 10
void IsotropicMaterial::computeDPFhat2DFij(const float *U, const float *V, const float * hessian, int i, int j, float *dPFhatdFij_diagonal)
{
	float w[3];
	// w[k] = dlambda_k/dF_ij = U_ik * V_jk see equation (7) of paper [PAPADOPOULO 2006]
	// "Estimating the Jacobian of the Singular Value Decomposition: Theory and Applications" 
	for (int k = 0; k < 3; ++k)
	{
		w[k] = U[i*3 + k] * V[j*3 + k];
	}

	for (int k = 0; k < 3; ++k)
	{
		dPFhatdFij_diagonal[k] = hessian[k] * w[0] + hessian[k + 3] * w[1] + hessian[k + 6] * w[2];
	}
}

void IsotropicMaterial::computeDP2DF(int tetID, const float *U, const float *Fhats, const float *V, float *dPdF)
{
	Eigen::Matrix3f Utilde, Vtilde;
	Utilde << U[0], U[1], U[2],
		U[3], U[4], U[5],
		U[6], U[7], U[8];
	Vtilde << V[0], V[1], V[2],
		V[3], V[4], V[5],
		V[6], V[7], V[8];

	float Ftildes[3];

	// perturbation: handle degenerated cases: see the paragraph between equation (9) and equation (10)
	// in paper [Xu et al. 2015]
	bool isPerturbed = true;
	float eps_singularvalue = 1e-6;
	// attention: Fhats are already sorted in descending order.
	if (Fhats[0] - Fhats[2] < 2 * eps_singularvalue)
	{
		Ftildes[2] = Fhats[2];
		Ftildes[1] = Ftildes[2] + eps_singularvalue;
		Ftildes[0] = Ftildes[1] + eps_singularvalue;
	}
	else // Fhats[0] - Fhats[2] >= 2 * m_eps_singularvalue
	{
		if ((Fhats[0] - Fhats[1] < eps_singularvalue) && (Fhats[1] - Fhats[2] >= eps_singularvalue))
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[0] - eps_singularvalue;
			Ftildes[0] = Fhats[0];
		}
		else if ((Fhats[0] - Fhats[1] >= eps_singularvalue) && (Fhats[1] - Fhats[2] < eps_singularvalue))
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[2] + eps_singularvalue;
			Ftildes[0] = Fhats[0];
		}
		else
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[1];
			Ftildes[0] = Fhats[0];
			isPerturbed = false;
		}
	}

	Eigen::Matrix3f Fnew, Unew, Vnew;
	Eigen::Vector3f  Fhatnew;

	if (isPerturbed)
	{
		Utilde.col(0) *= Ftildes[0];
		Utilde.col(1) *= Ftildes[1];
		Utilde.col(2) *= Ftildes[2];
		Fnew = Utilde * Vtilde.transpose();
		computeSVD33modified(Fnew, Fhatnew, Unew, Vnew);
	}
	else
	{
		for (int i = 0; i < 3;++i)
			Fhatnew[i] = Fhats[i];
		Unew = Utilde;
		Vnew = Vtilde;
	}
	//compute PFhat
	float PFhat[3];
	computeEnergy2FhatGradient(tetID, Fhatnew.data(), PFhat);

	float hessian[9];
	computeEnergy2FhatHessian(tetID, Fhatnew.data(), hessian);


	float dPdFij[9];
	int Fid;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			Fid = 3 * i + j;
			// transpose Unew and Vnew cause they are in stored in column-major matrix (Eigen default setting)
			computeDPDFij(Unew.transpose().data(), Fhatnew.data(), Vnew.transpose().data(),
				           PFhat, hessian, i, j, dPdFij);
			// copy dPdFij to dPdF
			for (int k = 0; k < 9; ++k)
			{
				dPdF[k * 9 + Fid] = dPdFij[k];
			}
		}
	}
}

void IsotropicMaterial::computeDPDFij(const float *U, const float *Fhats, const float *V, const float *PFhats, const float *hessian, int i, int j, float *dPdFij)
{
	Eigen::Matrix3f Umat, Vmat;
	Umat << U[0], U[1], U[2],
		U[3], U[4], U[5],
		U[6], U[7], U[8];
	Vmat << V[0], V[1], V[2],
		V[3], V[4], V[5],
		V[6], V[7], V[8];

	Eigen::Matrix3f wU, wVT;
	wU = Eigen::Matrix3f::Zero();
	wVT = Eigen::Matrix3f::Zero();

	Eigen::Matrix2d lambdaMat;
	Eigen::Vector2d uv,wUVT;

	int kl[3][2] = { { 0, 1 }, { 0, 2 }, { 1, 2 } };
	int k, l;

	for (int kli = 0; kli < 3; ++kli)
	{
		k = kl[kli][0];
		l = kl[kli][1];
		lambdaMat(0, 0) = lambdaMat(1, 1) = Fhats[l];
		lambdaMat(0, 1) = lambdaMat(1, 0) = Fhats[k];
		uv[0] = U[3 * i + k] * V[3 * j + l];
		uv[1] = -U[3 * i + l] * V[3 * j + k];

		wUVT = lambdaMat.inverse() * uv;

		wU(k, l) = wUVT(0);
		wU(l, k) = -wUVT(0);

		wVT(k, l) = wUVT(1);
		wVT(l, k) = -wUVT(1);
	}

	Eigen::Matrix3f dUdFij = Umat * wU;
	Eigen::Matrix3f dVTdFij = wVT * Vmat.transpose();

	float dPFhatdFij[3];

	computeDPFhat2DFij(U, V, hessian, i, j, dPFhatdFij);

	// equation(5) in paper [Xu et al. 2015] section 3.2
	Eigen::Matrix3f dPdFijMat = helperMatDiagonalMat(dUdFij, PFhats, Vmat.transpose())
		+ helperMatDiagonalMat(Umat, dPFhatdFij, Vmat.transpose())
		+ helperMatDiagonalMat(Umat, PFhats, dVTdFij);

	for (int m = 0; m < 3; ++m)
	{
		for (int n = 0; n < 3; ++n)
		{
			dPdFij[m * 3 + n] = dPdFijMat(m, n);
		}
	}
}

Eigen::Matrix3f IsotropicMaterial::helperMatDiagonalMat(Eigen::Matrix3f A, const float *diagonal, Eigen::Matrix3f B)
{
	for (int i = 0; i < 3;++i)
		A.col(i) *= diagonal[i];

	return A*B;
}