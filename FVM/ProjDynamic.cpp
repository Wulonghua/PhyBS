#include "ProjDynamic.h"


ProjDynamic::ProjDynamic(std::shared_ptr<TetMesh> tetMesh)
{
	n_nodes = tetMesh->getNodesNum();
	n_tets = tetMesh->getTetsNum();
	m_globalSolverMat.resize(3 * n_nodes, 3 * n_nodes);
	m_stiffWeight.resize(n_tets);
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
		wAcTAc = m_stiffWeight * (Ac.transpose()*Ac);
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
}
