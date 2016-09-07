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

	for (size_t i = 0; i < n_tets; ++i)
	{
		Eigen::Matrix3d Dm = Eigen::Matrix3d::Zero();
		Dm.col(0) = (m_nodes.col(m_tets(1, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(1) = (m_nodes.col(m_tets(2, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(2) = (m_nodes.col(m_tets(3, i)) - m_nodes.col(m_tets(0, i)));
		m_Dm_inverses.block<3,3>(0,i*3) = Dm.inverse();
	}

	std::cout << "tet model has been initialized.";
}

void TetMesh::SetTetMaterial(double e, double nu)
{
	m_E  = e;
	m_nu = nu;
}