#include "ProjDynamic.h"


ProjDynamic::ProjDynamic(std::shared_ptr<TetMesh> tetMesh):
m_stiffness(1e6), m_iterations(5)
{
	n_nodes = tetMesh->getNodesNum();
	n_tets = tetMesh->getTetsNum();
	//m_pos_new = Eigen::MatrixXf::Zero(3, n_nodes);
	m_pos_new = tetMesh->getNodes();
	m_globalSolverMat.resize(3 * n_nodes, 3 * n_nodes);
	m_stiffWeight.resize(n_tets);

	for (int i = 0; i < n_tets; ++i)
		m_stiffWeight[i] = m_stiffness * tetMesh->getTetVolumes()[i];

	initLocalProjection(tetMesh->getTets());
}


ProjDynamic::~ProjDynamic()
{
}

void ProjDynamic::initLocalProjection(const Eigen::MatrixXi & tets)
{

	m_localProjections = Eigen::MatrixXf::Zero(12, n_tets);
	m_proj_idx.resize(n_nodes);
	
	for (int i = 0; i < n_nodes; ++i)
	{
		m_proj_idx[i].reserve(4);
	}

	for (int i = 0; i < n_tets; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			m_proj_idx[tets(j, i)].push_back(i * 12 + j * 3);
		}
	}

	std::cout << "Projective Dynamics initialized." << std::endl;
}


void ProjDynamic::buildGlobalSolverMatrix(const Eigen::VectorXf &node_mass, const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Dm_inverse)
{
	typedef Eigen::Triplet<float> Tri;
	std::vector<Tri> triList;
	//triList.reserve(3 * n_nodes + 144*n_tets);
	float t2_inv = 1 / (t *t);
	for (int i = 0; i < n_nodes; ++i)
	{
		int ii = i * 3;
		triList.push_back(Tri(ii,   ii,   node_mass[i] * t2_inv));
		triList.push_back(Tri(ii+1, ii+1, node_mass[i] * t2_inv));
		triList.push_back(Tri(ii+2, ii+2, node_mass[i] * t2_inv));
	}

	Eigen::Matrix3f DmInvT;
	Eigen::Vector3f tmp;
	Eigen::MatrixXf Ac(9, 12);
	Eigen::MatrixXf wAcTAc(12, 12);
	Eigen::Vector4i v;  // nodes' indices
	int lRow, lCol, gRow, gCol;
	for (int k = 0; k < n_tets; ++k)
	{
		DmInvT = Dm_inverse.block<3, 3>(0, k * 3).transpose();
		v = tets.col(k);
		tmp = -1.0 * DmInvT.rowwise().sum();
		Ac.setZero();
		Ac.block<3, 1>(0, 0) = Ac.block<3, 1>(3, 1) = Ac.block<3, 1>(6, 2) = tmp;
		Ac.block<3, 1>(0, 3) = Ac.block<3, 1>(3, 4) = Ac.block<3, 1>(6, 5) = DmInvT.col(0);
		Ac.block<3, 1>(0, 6) = Ac.block<3, 1>(3, 7) = Ac.block<3, 1>(6, 8) = DmInvT.col(1);
		Ac.block<3, 1>(0, 9) = Ac.block<3, 1>(3, 10) = Ac.block<3, 1>(6, 11) = DmInvT.col(2);
		wAcTAc = m_stiffWeight[k] * (Ac.transpose()*Ac);
		for (int i = 0; i < 4; ++i)
		{
			for (int ii = 0; ii < 3; ++ii)
			{
				lRow = i * 3 + ii;
				gRow = v[i] * 3 + ii;
				for (int j = 0; j < 4; ++j)
				{
					for (int jj = 0; jj < 3; ++jj)
					{
						lCol = j * 3 + jj;
						gCol = v[j] * 3 + jj;
						if (wAcTAc(lRow, lCol) == 0) continue;
						triList.push_back(Tri(gRow, gCol, wAcTAc(lRow,lCol)));
					}
				}
			}
		}
	}
	m_globalSolverMat.setFromTriplets(triList.begin(), triList.end());
	m_pardiso_solver.compute(m_globalSolverMat);
}

Eigen::VectorXf ProjDynamic::projectLocalConstraints(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
													 const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, 
													 const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
													 const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext)
{
	
	Eigen::Matrix3f Ds, F, U, R, V, DmInv, DmInvT;
	//Eigen::Vector3f tmp;
	float tmp0, tmp1, tmp2;
	Eigen::Vector4i v;  // nodes' indices
	Eigen::VectorXf Rv,c;
	Rv.resize(9);
	c.resize(12);

	for (int k = 0; k < n_tets; ++k)
	{

		Ds.col(0) = pos.col(tets(1, k)) - pos.col(tets(0, k));
		Ds.col(1) = pos.col(tets(2, k)) - pos.col(tets(0, k));
		Ds.col(2) = pos.col(tets(3, k)) - pos.col(tets(0, k));
		DmInv = Dm_inverse.block<3, 3>(0, 3 * k);
		F = Ds * DmInv;
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		if (U.determinant() < 0)
			U.col(0) = -U.col(0);
		if (V.determinant() < 0)
			V.col(0) = -V.col(0);

		R = U * V.transpose();


		Rv(0) = R(0, 0); Rv(1) = R(0, 1); Rv(2) = R(0, 2);
		Rv(3) = R(1, 0); Rv(4) = R(1, 1); Rv(5) = R(1, 2);
		Rv(6) = R(2, 0); Rv(7) = R(2, 1); Rv(8) = R(2, 2);

		v = tets.col(k);

		//DmInvT = Dm_inverse.block<3, 3>(0, k * 3).transpose();
		//tmp = -1.0 * DmInvT.rowwise().sum();
		//Eigen::MatrixXf AcT(9, 12);
		//AcT.setZero();
		//AcT.block<3, 1>(0, 0) = AcT.block<3, 1>(3, 1)  = AcT.block<3, 1>(6, 2)  = tmp;
		//AcT.block<3, 1>(0, 3) = AcT.block<3, 1>(3, 4)  = AcT.block<3, 1>(6, 5)  = DmInvT.col(0);
		//AcT.block<3, 1>(0, 6) = AcT.block<3, 1>(3, 7)  = AcT.block<3, 1>(6, 8)  = DmInvT.col(1);
		//AcT.block<3, 1>(0, 9) = AcT.block<3, 1>(3, 10) = AcT.block<3, 1>(6, 11) = DmInvT.col(2);
		//AcT.transposeInPlace(); // this function is rather slow
		
		//Eigen::MatrixXf AcT(12, 9);
		//AcT.setZero();
		tmp0 = -DmInv(0, 0) - DmInv(1, 0) - DmInv(2, 0);
		tmp1 = -DmInv(0, 1) - DmInv(1, 1) - DmInv(2, 1);
		tmp2 = -DmInv(0, 2) - DmInv(1, 2) - DmInv(2, 2);
		//AcT(0, 0) = AcT(1, 3) = AcT(2, 6) = tmp0;
		//AcT(0, 1) = AcT(1, 4) = AcT(2, 7) = tmp1;
		//AcT(0, 2) = AcT(1, 5) = AcT(2, 8) = tmp2;
		//AcT(3, 0) = AcT(4, 3) = AcT(5, 6) = DmInv(0, 0);
		//AcT(3, 1) = AcT(4, 4) = AcT(5, 7) = DmInv(0, 1);
		//AcT(3, 2) = AcT(4, 5) = AcT(5, 8) = DmInv(0, 2);
		//AcT(6, 0) = AcT(7, 3) = AcT(8, 6) = DmInv(1, 0);
		//AcT(6, 1) = AcT(7, 4) = AcT(8, 7) = DmInv(1, 1);
		//AcT(6, 2) = AcT(7, 5) = AcT(8, 8) = DmInv(1, 2);
		//AcT(9, 0) = AcT(10, 3) = AcT(11, 6) = DmInv(2, 0);
		//AcT(9, 1) = AcT(10, 4) = AcT(11, 7) = DmInv(2, 1);
		//AcT(9, 2) = AcT(10, 5) = AcT(11, 8) = DmInv(2, 2);

		//c = m_stiffWeight[k] * AcT * Rv;
		c[0] = tmp0 * Rv[0] + tmp1 * Rv[1] + tmp2 * Rv[2];
		c[1] = tmp0 * Rv[3] + tmp1 * Rv[4] + tmp2 * Rv[5];
		c[2] = tmp0 * Rv[6] + tmp1 * Rv[7] + tmp2 * Rv[8];
		c[3] = DmInv(0, 0) * Rv[0] + DmInv(0, 1) * Rv[1] + DmInv(0, 2) * Rv[2];
		c[4] = DmInv(0, 0) * Rv[3] + DmInv(0, 1) * Rv[4] + DmInv(0, 2) * Rv[5];
		c[5] = DmInv(0, 0) * Rv[6] + DmInv(0, 1) * Rv[7] + DmInv(0, 2) * Rv[8];
		c[6] = DmInv(1, 0) * Rv[0] + DmInv(1, 1) * Rv[1] + DmInv(1, 2) * Rv[2];
		c[7] = DmInv(1, 0) * Rv[3] + DmInv(1, 1) * Rv[4] + DmInv(1, 2) * Rv[5];
		c[8] = DmInv(1, 0) * Rv[6] + DmInv(1, 1) * Rv[7] + DmInv(1, 2) * Rv[8];
		c[9] = DmInv(2, 0) * Rv[0] + DmInv(2, 1) * Rv[1] + DmInv(2, 2) * Rv[2];
		c[10] = DmInv(2, 0) * Rv[3] + DmInv(2, 1) * Rv[4] + DmInv(2, 2) * Rv[5];
		c[11] = DmInv(2, 0) * Rv[6] + DmInv(2, 1) * Rv[7] + DmInv(2, 2) * Rv[8];

		c = c*m_stiffWeight[k];
		
		for (int i = 0; i < 4; ++i)
		{
			s.col(v(i)) += c.block<3, 1>(i * 3, 0); 
		}

	}

	s.resize(3*n_nodes,1);
	return s;
}

Eigen::VectorXf ProjDynamic::projectLocalConstraints2(const Eigen::VectorXf &node_mass, const Eigen::VectorXf &node_inv_mass,
													  const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s,
													  const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
													  const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext)
{

	Eigen::Matrix3f Ds, F, U, R, V, DmInv, DmInvT;
	//Eigen::Vector3f tmp;
	float tmp0, tmp1, tmp2;
	Eigen::Vector4i v;  // nodes' indices
	Eigen::VectorXf Rv, c;
	Rv.resize(9);
	c.resize(12);

	for (int k = 0; k < n_tets; ++k)
	{

		Ds.col(0) = pos.col(tets(1, k)) - pos.col(tets(0, k));
		Ds.col(1) = pos.col(tets(2, k)) - pos.col(tets(0, k));
		Ds.col(2) = pos.col(tets(3, k)) - pos.col(tets(0, k));
		DmInv = Dm_inverse.block<3, 3>(0, 3 * k);
		F = Ds * DmInv;
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		if (U.determinant() < 0)
			U.col(0) = -U.col(0);
		if (V.determinant() < 0)
			V.col(0) = -V.col(0);

		R = U * V.transpose();


		Rv(0) = R(0, 0); Rv(1) = R(0, 1); Rv(2) = R(0, 2);
		Rv(3) = R(1, 0); Rv(4) = R(1, 1); Rv(5) = R(1, 2);
		Rv(6) = R(2, 0); Rv(7) = R(2, 1); Rv(8) = R(2, 2);

		v = tets.col(k);

		tmp0 = -DmInv(0, 0) - DmInv(1, 0) - DmInv(2, 0);
		tmp1 = -DmInv(0, 1) - DmInv(1, 1) - DmInv(2, 1);
		tmp2 = -DmInv(0, 2) - DmInv(1, 2) - DmInv(2, 2);

		c[0] = tmp0 * Rv[0] + tmp1 * Rv[1] + tmp2 * Rv[2];
		c[1] = tmp0 * Rv[3] + tmp1 * Rv[4] + tmp2 * Rv[5];
		c[2] = tmp0 * Rv[6] + tmp1 * Rv[7] + tmp2 * Rv[8];
		c[3] = DmInv(0, 0) * Rv[0] + DmInv(0, 1) * Rv[1] + DmInv(0, 2) * Rv[2];
		c[4] = DmInv(0, 0) * Rv[3] + DmInv(0, 1) * Rv[4] + DmInv(0, 2) * Rv[5];
		c[5] = DmInv(0, 0) * Rv[6] + DmInv(0, 1) * Rv[7] + DmInv(0, 2) * Rv[8];
		c[6] = DmInv(1, 0) * Rv[0] + DmInv(1, 1) * Rv[1] + DmInv(1, 2) * Rv[2];
		c[7] = DmInv(1, 0) * Rv[3] + DmInv(1, 1) * Rv[4] + DmInv(1, 2) * Rv[5];
		c[8] = DmInv(1, 0) * Rv[6] + DmInv(1, 1) * Rv[7] + DmInv(1, 2) * Rv[8];
		c[9] = DmInv(2, 0) * Rv[0] + DmInv(2, 1) * Rv[1] + DmInv(2, 2) * Rv[2];
		c[10] = DmInv(2, 0) * Rv[3] + DmInv(2, 1) * Rv[4] + DmInv(2, 2) * Rv[5];
		c[11] = DmInv(2, 0) * Rv[6] + DmInv(2, 1) * Rv[7] + DmInv(2, 2) * Rv[8];

		c = c*m_stiffWeight[k];
		
		m_localProjections.col(k) = c;
	}

	//for (int i = 0; i < 4; ++i)
	//{
	//	s.col(v(i)) += c.block<3, 1>(i * 3, 0);
	//}
	Eigen::Map<Eigen::VectorXf> proj(m_localProjections.data(), n_tets * 12);
	for (int i = 0; i < n_nodes; ++i)
	{
		int n = m_proj_idx[i].size();
		for (int j = 0; j < n; ++j)
		{
			s(0, i) += proj[m_proj_idx[i][j]];
			s(1, i) += proj[m_proj_idx[i][j] + 1];
			s(2, i) += proj[m_proj_idx[i][j] + 2];
		}
	}

	s.resize(3 * n_nodes, 1);
	return s;
}

Eigen::VectorXf ProjDynamic::projectLocalConstraints(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
	const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s,
	const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
	const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext, int num_threads)
{
	omp_set_num_threads(num_threads);
#pragma omp parallel for
	for (int k = 0; k < n_tets; ++k)
	{
		Eigen::Matrix3f Ds, F, U, R, V, DmInv;
		float tmp0, tmp1, tmp2;
		Eigen::Vector4i v;  // nodes' indices
		Eigen::VectorXf Rv, c;
		Rv.resize(9);
		c.resize(12);

		Ds.col(0) = pos.col(tets(1, k)) - pos.col(tets(0, k));
		Ds.col(1) = pos.col(tets(2, k)) - pos.col(tets(0, k));
		Ds.col(2) = pos.col(tets(3, k)) - pos.col(tets(0, k));
		DmInv = Dm_inverse.block<3, 3>(0, 3 * k);
		F = Ds * DmInv;
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		if (U.determinant() < 0)
			U.col(0) = -U.col(0);
		if (V.determinant() < 0)
			V.col(0) = -V.col(0);

		R = U * V.transpose();


		Rv(0) = R(0, 0); Rv(1) = R(0, 1); Rv(2) = R(0, 2);
		Rv(3) = R(1, 0); Rv(4) = R(1, 1); Rv(5) = R(1, 2);
		Rv(6) = R(2, 0); Rv(7) = R(2, 1); Rv(8) = R(2, 2);

		v = tets.col(k);

		tmp0 = -DmInv(0, 0) - DmInv(1, 0) - DmInv(2, 0);
		tmp1 = -DmInv(0, 1) - DmInv(1, 1) - DmInv(2, 1);
		tmp2 = -DmInv(0, 2) - DmInv(1, 2) - DmInv(2, 2);

		c[0] = tmp0 * Rv[0] + tmp1 * Rv[1] + tmp2 * Rv[2];
		c[1] = tmp0 * Rv[3] + tmp1 * Rv[4] + tmp2 * Rv[5];
		c[2] = tmp0 * Rv[6] + tmp1 * Rv[7] + tmp2 * Rv[8];
		c[3] = DmInv(0, 0) * Rv[0] + DmInv(0, 1) * Rv[1] + DmInv(0, 2) * Rv[2];
		c[4] = DmInv(0, 0) * Rv[3] + DmInv(0, 1) * Rv[4] + DmInv(0, 2) * Rv[5];
		c[5] = DmInv(0, 0) * Rv[6] + DmInv(0, 1) * Rv[7] + DmInv(0, 2) * Rv[8];
		c[6] = DmInv(1, 0) * Rv[0] + DmInv(1, 1) * Rv[1] + DmInv(1, 2) * Rv[2];
		c[7] = DmInv(1, 0) * Rv[3] + DmInv(1, 1) * Rv[4] + DmInv(1, 2) * Rv[5];
		c[8] = DmInv(1, 0) * Rv[6] + DmInv(1, 1) * Rv[7] + DmInv(1, 2) * Rv[8];
		c[9] = DmInv(2, 0) * Rv[0] + DmInv(2, 1) * Rv[1] + DmInv(2, 2) * Rv[2];
		c[10] = DmInv(2, 0) * Rv[3] + DmInv(2, 1) * Rv[4] + DmInv(2, 2) * Rv[5];
		c[11] = DmInv(2, 0) * Rv[6] + DmInv(2, 1) * Rv[7] + DmInv(2, 2) * Rv[8];

		c = c*m_stiffWeight[k];

		for (int i = 0; i < 4; ++i)
		{
#pragma omp atomic
			s(0, v(i)) += c(i * 3);
#pragma omp atomic
			s(1, v(i)) += c(i * 3 + 1);
#pragma omp atomic
			s(2, v(i)) += c(i * 3 + 2);
			//s.col(v(i)) += c.block<3, 1>(i * 3, 0);
		}

	}

	s.resize(3 * n_nodes, 1);
	return s;
}

Eigen::VectorXf ProjDynamic::projectLocalConstraints2(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
	const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
	const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext, int num_threads)
{

	omp_set_num_threads(num_threads);
#pragma omp parallel for
	for (int k = 0; k < n_tets; ++k)
	{
		Eigen::Matrix3f Ds, F, U, R, V, DmInv;
		//Eigen::Vector3f tmp;
		float tmp0, tmp1, tmp2;
		Eigen::Vector4i v;  // nodes' indices
		Eigen::VectorXf Rv, c;
		Rv.resize(9);
		c.resize(12);

		Ds.col(0) = pos.col(tets(1, k)) - pos.col(tets(0, k));
		Ds.col(1) = pos.col(tets(2, k)) - pos.col(tets(0, k));
		Ds.col(2) = pos.col(tets(3, k)) - pos.col(tets(0, k));
		DmInv = Dm_inverse.block<3, 3>(0, 3 * k);
		F = Ds * DmInv;
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		U = svd.matrixU();
		V = svd.matrixV();
		if (U.determinant() < 0)
			U.col(0) = -U.col(0);
		if (V.determinant() < 0)
			V.col(0) = -V.col(0);

		R = U * V.transpose();


		Rv(0) = R(0, 0); Rv(1) = R(0, 1); Rv(2) = R(0, 2);
		Rv(3) = R(1, 0); Rv(4) = R(1, 1); Rv(5) = R(1, 2);
		Rv(6) = R(2, 0); Rv(7) = R(2, 1); Rv(8) = R(2, 2);

		v = tets.col(k);

		tmp0 = -DmInv(0, 0) - DmInv(1, 0) - DmInv(2, 0);
		tmp1 = -DmInv(0, 1) - DmInv(1, 1) - DmInv(2, 1);
		tmp2 = -DmInv(0, 2) - DmInv(1, 2) - DmInv(2, 2);

		c[0] = tmp0 * Rv[0] + tmp1 * Rv[1] + tmp2 * Rv[2];
		c[1] = tmp0 * Rv[3] + tmp1 * Rv[4] + tmp2 * Rv[5];
		c[2] = tmp0 * Rv[6] + tmp1 * Rv[7] + tmp2 * Rv[8];
		c[3] = DmInv(0, 0) * Rv[0] + DmInv(0, 1) * Rv[1] + DmInv(0, 2) * Rv[2];
		c[4] = DmInv(0, 0) * Rv[3] + DmInv(0, 1) * Rv[4] + DmInv(0, 2) * Rv[5];
		c[5] = DmInv(0, 0) * Rv[6] + DmInv(0, 1) * Rv[7] + DmInv(0, 2) * Rv[8];
		c[6] = DmInv(1, 0) * Rv[0] + DmInv(1, 1) * Rv[1] + DmInv(1, 2) * Rv[2];
		c[7] = DmInv(1, 0) * Rv[3] + DmInv(1, 1) * Rv[4] + DmInv(1, 2) * Rv[5];
		c[8] = DmInv(1, 0) * Rv[6] + DmInv(1, 1) * Rv[7] + DmInv(1, 2) * Rv[8];
		c[9] = DmInv(2, 0) * Rv[0] + DmInv(2, 1) * Rv[1] + DmInv(2, 2) * Rv[2];
		c[10] = DmInv(2, 0) * Rv[3] + DmInv(2, 1) * Rv[4] + DmInv(2, 2) * Rv[5];
		c[11] = DmInv(2, 0) * Rv[6] + DmInv(2, 1) * Rv[7] + DmInv(2, 2) * Rv[8];

		c = c*m_stiffWeight[k];

		m_localProjections.col(k) = c;
	}

	//for (int i = 0; i < 4; ++i)
	//{
	//	s.col(v(i)) += c.block<3, 1>(i * 3, 0);
	//}
	Eigen::Map<Eigen::VectorXf> proj(m_localProjections.data(), n_tets * 12);
	
	omp_set_num_threads(num_threads);
#pragma omp parallel for
	for (int i = 0; i < n_nodes; ++i)
	{
		int n = m_proj_idx[i].size();
		for (int j = 0; j < n; ++j)
		{
			s(0, i) += proj[m_proj_idx[i][j]];
			s(1, i) += proj[m_proj_idx[i][j] + 1];
			s(2, i) += proj[m_proj_idx[i][j] + 2];
		}
	}

	s.resize(3 * n_nodes, 1);
	return s;
}


void ProjDynamic::solveGlobalStep(Eigen::MatrixXf &pos, Eigen::VectorXf &b)
{
	Eigen::Map<Eigen::VectorXf> p(m_pos_new.data(), n_nodes * 3);
	p = m_pardiso_solver.solve(b);
}

void ProjDynamic::doProjDynamics(Eigen::MatrixXf &pos, Eigen::MatrixXf &vel,
	const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
	const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Dm_inverse, const Eigen::MatrixXf & fext)
{
	m_pos_new = pos;
	Eigen::VectorXf b;

	Eigen::MatrixXf s = Eigen::MatrixXf::Zero(3, n_nodes);
	for (int i = 0; i < n_nodes; ++i)
	{
		s.col(i) = t*t * node_inv_mass[i] * fext.col(i);
	}

	s = pos + t * vel + s;

	for (int i = 0; i < n_nodes; ++i)
	{
		s.col(i) = s.col(i)*node_mass[i] / (t*t);
	}

	//m_time.start();
	//m_elapses1 = m_elapses2 = 0;
	for (int i = 0; i < m_iterations; ++i)
	{
		//m_time.restart();
		//b = projectLocalConstraints(node_mass, node_inv_mass, tets, t, s,m_pos_new, Dm_inverse, vel, fext,4);

		b = projectLocalConstraints2(node_mass, node_inv_mass, tets, t, s, m_pos_new, Dm_inverse, vel, fext, 4);
		//m_elapses1 += m_time.restart();
		solveGlobalStep(m_pos_new, b);
		//m_elapses2 += m_time.restart();
	}
	//m_elapses1 += m_time.restart();
	//std::cout << "time : " << m_elapses1 << std::endl;
	//std::cout << "time for local AcT:" << m_elapses2 << std::endl;
	m_time.invalidate();
	vel = (m_pos_new - pos) / t *0.99;
	pos = m_pos_new;
}
