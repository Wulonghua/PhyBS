#include "ProjDynamic.h"


ProjDynamic::ProjDynamic(std::shared_ptr<TetMesh> tetMesh):
m_stiffness(5000), m_iterations(50)
{
	n_nodes = tetMesh->getNodesNum();
	n_tets = tetMesh->getTetsNum();
	//m_pos_new = Eigen::MatrixXf::Zero(3, n_nodes);
	m_pos_new = tetMesh->getNodes();

	m_globalSolverMat.resize(3 * n_nodes, 3 * n_nodes);
	m_stiffWeight.resize(n_tets);

	for (int i = 0; i < n_tets; ++i)
		m_stiffWeight[i] = m_stiffness * tetMesh->getTetVolumes()[i];
}


ProjDynamic::~ProjDynamic()
{
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
													 const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
													 const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext)
{
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
	
	Eigen::Matrix3f Ds, F, DmInvT, U, R, V;
	Eigen::Vector3f tmp;
	Eigen::Vector4i v;  // nodes' indices
	Eigen::VectorXf Rv,c;
	Rv.resize(9);

	for (int k = 0; k < n_tets; ++k)
	{
		Ds.col(0) = pos.col(tets(1, k)) - pos.col(tets(0, k));
		Ds.col(1) = pos.col(tets(2, k)) - pos.col(tets(0, k));
		Ds.col(2) = pos.col(tets(3, k)) - pos.col(tets(0, k));

		F = Ds * Dm_inverse.block<3, 3>(0, 3*k);
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

		DmInvT = Dm_inverse.block<3, 3>(0, k * 3).transpose();
		v = tets.col(k);
		tmp = -1.0 * DmInvT.rowwise().sum();
		Eigen::MatrixXf AcT(9, 12);
		AcT.setZero();
		AcT.block<3, 1>(0, 0) = AcT.block<3, 1>(3, 1)  = AcT.block<3, 1>(6, 2)  = tmp;
		AcT.block<3, 1>(0, 3) = AcT.block<3, 1>(3, 4)  = AcT.block<3, 1>(6, 5)  = DmInvT.col(0);
		AcT.block<3, 1>(0, 6) = AcT.block<3, 1>(3, 7)  = AcT.block<3, 1>(6, 8)  = DmInvT.col(1);
		AcT.block<3, 1>(0, 9) = AcT.block<3, 1>(3, 10) = AcT.block<3, 1>(6, 11) = DmInvT.col(2);
		AcT.transposeInPlace();

		c = m_stiffWeight[k] * AcT * Rv;

		for (int i = 0; i < 4; ++i)
		{
			s.col(v(i)) += c.block<3, 1>(i * 3, 0);
		}
	}

	s.resize(3*n_nodes,1);
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
	for (int i = 0; i < m_iterations; ++i)
	{
		b = projectLocalConstraints(node_mass, node_inv_mass, tets, t, m_pos_new, Dm_inverse, vel, fext);
		solveGlobalStep(m_pos_new, b);
	}
	vel = (m_pos_new - pos) / t;
	pos = m_pos_new;
}
