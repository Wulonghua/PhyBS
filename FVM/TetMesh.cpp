#include "TetMesh.h"


TetMesh::TetMesh()
{
	InitModel();
}


TetMesh::~TetMesh()
{
}

void TetMesh::LoadNodesFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "Cannot load nodes file." << std::endl;
		exit(0);
	}
		
	QTextStream fin(&file);
	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_nodes = segs[0].toInt();
		
		//load nodes' positions
		int prefix;
		m_nodes = Eigen::MatrixXd::Zero(3, n_nodes);
		for (size_t i = 0; i < n_nodes; ++i)
		{
			fin >> prefix >> m_nodes(0, i) >> m_nodes(1, i) >> m_nodes(2, i);
		}
		//std::cout << m_nodes;
	}
}

void TetMesh::LoadTetsFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "Cannot load nodes file." << std::endl;
		exit(0);
	}
	QTextStream fin(&file);
	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_tets = segs[0].toInt();

		//load Tets' indices
		int prefix;
		m_tets = Eigen::MatrixXi::Zero(4, n_tets);
		for (size_t i = 0; i < n_tets; ++i)
		{
			fin >> prefix >> m_tets(0, i) >> m_tets(1, i) >> m_tets(2, i) >> m_tets(3, i);
		}
		//std::cout << m_tets;
	}
}

void TetMesh::LoadFacesFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cout << "Cannot load nodes file." << std::endl;
		exit(0);
	}
	QTextStream fin(&file);

	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_bound_faces = segs[0].toInt();

		//load Tets' indices
		int prefix;
		m_bound_faces = Eigen::MatrixXi::Zero(3, n_bound_faces);
		for (size_t i = 0; i < n_tets; ++i)
		{
			fin >> prefix >> m_bound_faces(0,i) >> m_bound_faces(1,i) >> m_bound_faces(2,i);
		}
		//std::cout << m_bound_faces;
	}
}

void TetMesh::InitModel()
{
	LoadNodesFromFile(QStringLiteral("..\\model\\box\\box.1.node"));
	LoadTetsFromFile(QStringLiteral("..\\model\\box\\box.1.ele"));
	LoadFacesFromFile(QStringLiteral("..\\model\\box\\box.1.face"));

	m_Dm_inverses = Eigen::MatrixXd::Zero(3, n_tets * 3);
	m_ANs		  = Eigen::MatrixXd::Zero(3, n_tets * 3);

	// precompute dim_inverse for each tetrahedron and bi for three nodes in each tetrahedron
	for (size_t i = 0; i < n_tets; ++i)
	{
		Eigen::Matrix3d Dm = Eigen::Matrix3d::Zero();
		Dm.col(0) = (m_nodes.col(m_tets(1, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(1) = (m_nodes.col(m_tets(2, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(2) = (m_nodes.col(m_tets(3, i)) - m_nodes.col(m_tets(0, i)));
		m_Dm_inverses.block<3,3>(0,i*3) = Dm.inverse();

		ComputeANs(i);
	}	

	std::cout << "tet model has been initialized."<<std::endl;
}

void TetMesh::ComputeANs(int tetid)
{
	int edge[3][4] = {{1,0,2,3},{2,0,1,3},{3,0,1,2}};
	Eigen::Vector3d v1, v2, v3;
	for (size_t i = 0; i < 3; ++i)
	{
		v1 = m_nodes.col(m_tets(edge[i][1], tetid)) - m_nodes.col(m_tets(edge[i][0], tetid));
		v2 = m_nodes.col(m_tets(edge[i][2], tetid)) - m_nodes.col(m_tets(edge[i][0], tetid));
		v3 = m_nodes.col(m_tets(edge[i][3], tetid)) - m_nodes.col(m_tets(edge[i][0], tetid));

		// Explanation here: the following line computes b_i = (-1/3)(A_1*N_1+A_2*N_2+A_3*N_3)
		// Note that cross product of two edges is twice the area of a face times the normal,
		// So we can simply add one sixth of -sigma times the cross product to each of the three nodes.
		m_ANs.col(tetid * 3 + i) = -(v1.cross(v2) + v2.cross(v3) + v3.cross(v1)) / 6.0;
		std::cout << m_ANs.col(tetid*3+i)<<std::endl;
	}
}

void TetMesh::SetTetMaterial(double e, double nu)
{
	m_E  = e;
	m_nu = nu;
}