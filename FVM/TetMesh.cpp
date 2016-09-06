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
		return;
	}
		
	QTextStream fin(&file);
	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_nodes = segs[0].toInt();
		
		//load nodes' positions
		int prefix;
		m_nodes = Eigen::MatrixXd::Zero(n_nodes,3);
		for (size_t i = 0; i < n_nodes; ++i)
		{
			fin >> prefix >> m_nodes(i, 0) >> m_nodes(i, 1) >> m_nodes(i, 2);
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
		return;
	}
	QTextStream fin(&file);
	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_tets = segs[0].toInt();

		//load Tets' indices
		int prefix;
		m_tets = Eigen::MatrixXi::Zero(n_tets, 4);
		for (size_t i = 0; i < n_tets; ++i)
		{
			fin >> prefix >> m_tets(i, 0) >> m_tets(i, 1) >> m_tets(i, 2) >> m_tets(i, 3);
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
		return;
	}
	QTextStream fin(&file);

	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_bound_faces = segs[0].toInt();

		//load Tets' indices
		int prefix;
		m_bound_faces = Eigen::MatrixXi::Zero(n_bound_faces, 3);
		for (size_t i = 0; i < n_tets; ++i)
		{
			fin >> prefix >> m_bound_faces(i, 0) >> m_bound_faces(i, 1) >> m_bound_faces(i, 2);
		}
		//std::cout << m_bound_faces;
	}
}

void TetMesh::InitModel()
{
	LoadNodesFromFile(QStringLiteral("G:\\wlh_code\\PhyBS\\model\\box\\box.1.node"));
	LoadTetsFromFile(QStringLiteral("G:\\wlh_code\\PhyBS\\model\\box\\box.1.ele"));
	LoadFacesFromFile(QStringLiteral("G:\\wlh_code\\PhyBS\\model\\box\\box.1.face"));

	m_Dm_inverses = Eigen::MatrixXd::Zero(n_tets*3,3);
	m_ANs = Eigen::MatrixXd::Zero(n_tets*3,3);

	for (size_t i = 0; i < n_tets; ++i)
	{
		Eigen::Matrix3d Dm = Eigen::Matrix3d::Zero();
		Dm.col(0) = (m_nodes.row(m_tets(i, 1)) - m_nodes.row(m_tets(i, 0))).transpose();
		Dm.col(1) = (m_nodes.row(m_tets(i, 2)) - m_nodes.row(m_tets(i, 0))).transpose();
		Dm.col(2) = (m_nodes.row(m_tets(i, 3)) - m_nodes.row(m_tets(i, 0))).transpose();
		m_Dm_inverses.block<3,3>(i*3,0) = Dm.inverse();
	}

	std::cout << "tet model has been initialized.";
}
