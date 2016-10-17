#include "IsotropicMaterial.h"


IsotropicMaterial::IsotropicMaterial() :
m_eps_singularvalue(1e-6), 
m_matrix33fromTeran({ { 0, 3, 5, 4, 1, 7, 6, 8, 2 } }) // compromised way to initialize: VS2013 does not fully support c++11
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
	//for (int i = 0; i < 3; ++i)
	//	S(i) = m_tetModel->fixPrecision(S(i));
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
	//Fix U if certain singularvalue is below epsilon or equal to 0.
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

	//Fix U and negate minimal singularvalue if determinant of U is equal to -1
	if (U.determinant() < 0) // ==-1
	{
		U.col(2) = U.col(2) * (-1);
		S(2) *= -1;
	}
}

void IsotropicMaterial::compouteFhats()
{
	int n = m_tetModel->getTetsNum();
	Eigen::Matrix3d F, U, V;
	Eigen::Vector3d Fhat;

	for (int i = 0; i < n; ++i)
	{
		F = m_tetModel->computeDeformationGradient(i);

		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;
		double t = Fhat[0];

		m_Us.block<3, 3>(0, i * 3) = U;
		m_Vs.block<3, 3>(0, i * 3) = V;
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

		//std::cout << "F: " << std::endl;
		//std::cout << F << std::endl;

		computeSVD33modified(F, Fhat, U, V);
		m_Fhats.col(i) = Fhat;
		double t = Fhat[0];

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

Eigen::MatrixXd IsotropicMaterial::computeDP2DF(int tetID)
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
					dPdFij_middle += UTeijV(i, j) * subTensor;
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

Eigen::MatrixXd IsotropicMaterial::computeInnerForcesfromFhats()
{
	computeFhatsInvariants();
	int n = m_tetModel->getTetsNum();
	Eigen::Vector3d Phat, Fhat, Fhat_inverseT;
	double mu, lambda, I3;
	Eigen::Matrix3d P, U, V, forces;

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



Eigen::MatrixXd IsotropicMaterial::computeStiffnessMatrix(int tetID)
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

	int convert_idx[9] = { 0, 3, 6, 1, 4, 7, 2, 5, 8 };
	for (int i = 0; i < 9; ++i)
	{
		dfdx.row(i + 3) = dGdx.row(convert_idx[i]);
	}

	// test
	//m_tetModel->writeMatrix("mat.csv", dfdx);
	return dfdx;
}

Eigen::SparseMatrix<double> IsotropicMaterial::computeGlobalStiffnessMatrix()
{
	int n = m_tetModel->getNodesNum();
	int m = m_tetModel->getTetsNum();
	Eigen::SparseMatrix<double> gK(3 * n, 3 * n);
	Eigen::MatrixXd K;
	int Ki, Kj, gKi, gKj;

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
						gK.coeffRef(gKi, gKj) += K(Ki, Kj);
					}
				}
			}
		}
	}

	return gK;
}

Eigen::Matrix3d IsotropicMaterial::restoreMatrix33fromTeranVector(Eigen::VectorXd v)
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

void IsotropicMaterial::computeEnergy2FhatGradient(int tetID, const double *Fhats, double *gradient)
{
	double lame_mu = m_mus[tetID];
	double lame_lambda = m_lambdas[tetID];
	double lambda1 = Fhats[0];
	double lambda2 = Fhats[1];
	double lambda3 = Fhats[2];

	double dg12 = dgEnergy(lambda1*lambda2, lame_mu, lame_lambda);
	double dg23 = dgEnergy(lambda2*lambda3, lame_mu, lame_lambda);
	double dg31 = dgEnergy(lambda3*lambda1, lame_mu, lame_lambda);

	double dh123 = dhEnergy(lambda1*lambda2*lambda3, lame_mu, lame_lambda);

	gradient[0] = dfEnergy(lambda1, lame_mu, lame_lambda) + dg12 * lambda2 + dg31 * lambda3 + dh123*lambda2*lambda3;
	gradient[1] = dfEnergy(lambda2, lame_mu, lame_lambda) + dg23 * lambda3 + dg12 * lambda1 + dh123*lambda3*lambda1;
	gradient[2] = dfEnergy(lambda3, lame_mu, lame_lambda) + dg31 * lambda1 + dg23 * lambda2 + dh123*lambda1*lambda2;
}


// *hessian is the one dimentional representation of the row-major 3*3 hessian matrix
void IsotropicMaterial::computeEnergy2FhatHessian(int tetID, const double *Fhats, double *hessian)
{
	double lame_mu = m_mus[tetID];
	double lame_lambda = m_lambdas[tetID];

	double lambda12 = Fhats[0] * Fhats[1];
	double lambda23 = Fhats[1] * Fhats[2];
	double lambda31 = Fhats[2] * Fhats[0];
	double lambda123 = lambda12 * Fhats[2];

	double dg12 = dgEnergy(lambda12, lame_mu, lame_lambda);
	double dg23 = dgEnergy(lambda23, lame_mu, lame_lambda);
	double dg31 = dgEnergy(lambda31, lame_mu, lame_lambda);

	double ddg12 = ddgEnergy(lambda12, lame_mu, lame_lambda);
	double ddg23 = ddgEnergy(lambda23, lame_mu, lame_lambda);
	double ddg31 = ddgEnergy(lambda31, lame_mu, lame_lambda);
	double dh123 = dhEnergy(lambda123, lame_mu, lame_lambda);
	double ddh123 = ddhEnergy(lambda123, lame_mu, lame_lambda);

	// hessian(1,1)
	hessian[0] = ddfEnergy(Fhats[0], lame_mu, lame_lambda) + ddg12 * Fhats[1] * Fhats[1]
														   + ddg31 * Fhats[2] * Fhats[2]
														   + ddh123 * lambda23 * lambda23;
	// hessian(2,2)
	hessian[4] = ddfEnergy(Fhats[1], lame_mu, lame_lambda) + ddg23 * Fhats[2] * Fhats[2]
														   + ddg12 * Fhats[0] * Fhats[0]
														   + ddh123 * lambda31 * lambda31;
	// hessian(3,3)
	hessian[8] = ddfEnergy(Fhats[2], lame_mu, lame_lambda) + ddg31 * Fhats[0] * Fhats[0]
														   + ddg23 * Fhats[1] * Fhats[1]
														   + ddh123 * lambda12 * lambda12;
	// hessian(1,2) = hessian(2,1)
	hessian[1] = hessian[3] = ddg12 * lambda12 + dg12 + ddh123 * lambda23 * lambda31 + dh123 * Fhats[2];

	// hessian(1,3) = hessian(3,1)
	hessian[2] = hessian[6] = ddg31 * lambda31 + dg31 + ddh123 * lambda12 * lambda23 + dh123 * Fhats[1];

	// hessian(2,3) = hessian(3,2)
	hessian[5] = hessian[7] = ddg23 * lambda23 + dg23 + ddh123 * lambda12 * lambda31 + dh123 * Fhats[0];
}

// see [Xu et al. 2015] Section 3.1 euqation 10
void IsotropicMaterial::computeDPFhat2DFij(const double *U, const double *V, const double * hessian, int i, int j, double *dPFhatdFij_diagonal)
{
	double w[3];
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

void IsotropicMaterial::computeDPDF(const double *U, const double *Fhats, const double *V)
{
	Eigen::Matrix3d Utilde, Vtilde;
	Utilde << U[0], U[1], U[2],
		U[3], U[4], U[5],
		U[6], U[7], U[8];
	Vtilde << V[0], V[1], V[2],
		V[3], V[4], V[5],
		V[6], V[7], V[8];

	double Ftildes[3];

	// perturbation: handle degenerated cases: see the paragraph between equation (9) and equation (10)
	// in paper [Xu et al. 2015]

	// attention: Fhats are already sorted in descending order.
	if (Fhats[0] - Fhats[2] < 2 * m_eps_singularvalue)
	{
		Ftildes[2] = Fhats[2];
		Ftildes[1] = Ftildes[2] + m_eps_singularvalue;
		Ftildes[0] = Ftildes[1] + m_eps_singularvalue;
	}
	else // Fhats[0] - Fhats[2] >= 2 * m_eps_singularvalue
	{
		if ((Fhats[0] - Fhats[1] < m_eps_singularvalue) && (Fhats[1] - Fhats[2] >= m_eps_singularvalue))
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[1];
			Ftildes[0] = Fhats[1] + m_eps_singularvalue;
		}
		else if ((Fhats[0] - Fhats[1] >= m_eps_singularvalue) && (Fhats[1] - Fhats[2] < m_eps_singularvalue))
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[2] + m_eps_singularvalue;
			Ftildes[0] = Fhats[0];
		}
		else
		{
			Ftildes[2] = Fhats[2];
			Ftildes[1] = Fhats[1];
			Ftildes[0] = Fhats[0];
		}
	}

	Utilde.col(0) *= Ftildes[0];
	Utilde.col(1) *= Ftildes[1];
	Utilde.col(2) *= Ftildes[2];

	Eigen::Matrix3d Fnew = Utilde * Vtilde.transpose();

	Eigen::Matrix3d Unew, Vnew;
	Eigen::Vector3d  Fhatnew;

	computeSVD33modified(Fnew, Fhatnew, Unew, Vnew);

}

void IsotropicMaterial::computeDPDFij(const double *U, const double *Fhat, const double *V)
{




}