#include "TetMesh.h"
#include "RenderWidget.h"


TetMesh::TetMesh()
{
	initModel();
}

TetMesh::~TetMesh()
{
}

void TetMesh::initNodesFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cerr << "Cannot load nodes file." << std::endl;
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
		m_rest_positions = Eigen::MatrixXd::Zero(3, n_nodes);

		for (size_t i = 0; i < n_nodes; ++i)
		{
			fin >> prefix >> m_rest_positions(0, i) >> m_rest_positions(1, i) >> m_rest_positions(2, i);

			//if (m_rest_positions(0, i) == -2.5)
			//	m_constraintIDs.push_back(i);

			//if (m_rest_positions(1, i) > 0.16)
			//	m_constraintIDs.push_back(i);

			if (m_rest_positions(1, i) < -1.8)
				m_constraintIDs.push_back(i);
		}

		m_nodes = m_rest_positions;


		m_nodes_mass = Eigen::VectorXd::Zero(n_nodes);
		m_velocities = Eigen::MatrixXd::Zero(3, n_nodes);

		//for test
		//m_velocities.row(1) = -50.0 * Eigen::VectorXd::Ones(n_nodes);

		m_nodes_gravity = Eigen::MatrixXd::Zero(3, n_nodes);
		m_nodes_forces = Eigen::MatrixXd::Zero(3, n_nodes);
		m_nodes_external_forces = Eigen::MatrixXd::Zero(3, n_nodes);
		//std::cout << m_nodes;
	}
}

void TetMesh::updateNodesFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cerr << "Cannot load nodes file." << std::endl;
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

		for (size_t i = 0; i < n_nodes; ++i)
		{
			fin >> prefix >> m_nodes(0, i) >> m_nodes(1, i) >> m_nodes(2, i);
		}
	}
}

void TetMesh::initTetsFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cerr << "Cannot load nodes file." << std::endl;
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

void TetMesh::initFacesFromFile(QString filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		std::cerr << "Cannot load nodes file." << std::endl;
		exit(1);
	}
	QTextStream fin(&file);

	if (!fin.atEnd())
	{
		QString line = fin.readLine().simplified();
		QStringList segs = line.split(" ");
		n_bound_faces = segs[0].toInt();
		m_bound_normals = Eigen::MatrixXd::Zero(3,n_bound_faces);
		m_face_centers = Eigen::MatrixXd::Zero(3,n_bound_faces);
		m_faceRingIndices.resize(n_bound_faces);

		//load faces' indices
		int prefix;
		m_bound_faces = Eigen::MatrixXi::Zero(3, n_bound_faces);
		for (size_t i = 0; i < n_bound_faces; ++i)
		{
			fin >> prefix >> m_bound_faces(0,i) >> m_bound_faces(1,i) >> m_bound_faces(2,i);
		}
	}
	computeBoundfaceRingIndices();
}

void TetMesh::initModel()
{
	//initNodesFromFile(QStringLiteral("..\\model\\torus\\torus.1.node"));
	//initTetsFromFile(QStringLiteral("..\\model\\torus\\torus.1.ele"));
	//initFacesFromFile(QStringLiteral("..\\model\\torus\\torus.1.face"));

	//initNodesFromFile(QStringLiteral("..\\model\\tet\\tet.1.node"));
	//initTetsFromFile(QStringLiteral("..\\model\\tet\\tet.1.ele"));
	//initFacesFromFile(QStringLiteral("..\\model\\tet\\tet.1.face"));

	//initNodesFromFile(QStringLiteral("..\\model\\bar\\bar.1.node"));
	//initTetsFromFile(QStringLiteral("..\\model\\bar\\bar.1.ele"));
	//initFacesFromFile(QStringLiteral("..\\model\\bar\\bar.1.face"));

	//initNodesFromFile(QStringLiteral("..\\model\\bar2\\bar2.1.node"));
	//initTetsFromFile(QStringLiteral("..\\model\\bar2\\bar2.1.ele"));
	//initFacesFromFile(QStringLiteral("..\\model\\bar2\\bar2.1.face"));

	//initNodesFromFile(QStringLiteral("..\\model\\bunny\\bunny.1.node"));
	//initTetsFromFile(QStringLiteral("..\\model\\bunny\\bunny.1.ele"));
	//initFacesFromFile(QStringLiteral("..\\model\\bunny\\bunny.1.face"));

	initNodesFromFile(QStringLiteral("..\\model\\asiandragon\\asiandragon.1.node"));
	initTetsFromFile(QStringLiteral("..\\model\\asiandragon\\asiandragon.1.ele"));
	initFacesFromFile(QStringLiteral("..\\model\\asiandragon\\asiandragon.1.face"));

	m_Dm_inverses = Eigen::MatrixXd::Zero(3, n_tets * 3);
	m_ANs		  = Eigen::MatrixXd::Zero(3, n_tets * 3);
	setTetMaterial(1000000, 0.45,1000);

	// precompute dim_inverse for each tetrahedron and bi for three nodes in each tetrahedron
	// also each node's weight
	double w;
	for (size_t i = 0; i < n_tets; ++i)
	{
		Eigen::Matrix3d Dm = Eigen::Matrix3d::Zero();
		Dm.col(0) = (m_nodes.col(m_tets(1, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(1) = (m_nodes.col(m_tets(2, i)) - m_nodes.col(m_tets(0, i)));
		Dm.col(2) = (m_nodes.col(m_tets(3, i)) - m_nodes.col(m_tets(0, i)));
		w = std::abs(Dm.determinant()) / 24.0 * m_density;
		for (size_t j = 0; j < 4; ++j)
		{
			m_nodes_mass(m_tets(j, i)) += w;
		}
		m_Dm_inverses.block<3,3>(0,i*3) = Dm.inverse();
		computeANs(i);
	}	
	m_nodes_gravity.row(1) = -9.8 * m_nodes_mass.transpose();
	
	std::cout << "tet model has been initialized."<<std::endl;
}

void TetMesh::computeBoundfaceRingIndices()
{
	std::vector<std::set<int>> nodeRing;
	std::set<int> tmp;
	int i1, i2, i3;

	nodeRing.resize(n_nodes);

	
	for (int i = 0; i < n_bound_faces; ++i)
	{
		nodeRing[m_bound_faces(0, i)].insert(m_bound_faces(1, i));
		nodeRing[m_bound_faces(0, i)].insert(m_bound_faces(2, i));
		nodeRing[m_bound_faces(1, i)].insert(m_bound_faces(0, i));
		nodeRing[m_bound_faces(1, i)].insert(m_bound_faces(2, i));
		nodeRing[m_bound_faces(2, i)].insert(m_bound_faces(0, i));
		nodeRing[m_bound_faces(2, i)].insert(m_bound_faces(1, i));
	}

	for (int i = 0; i < n_bound_faces; ++i)
	{
		tmp.clear();
		int ii[3];
		ii[0] = m_bound_faces(0, i);
		ii[1] = m_bound_faces(1, i);
		ii[2] = m_bound_faces(2, i);

		for (int j = 0; j < 3; ++j)
		{
			for (auto iter = nodeRing[ii[j]].begin(); iter != nodeRing[ii[j]].end(); ++iter)
			{
				tmp.insert(*iter);
			}
		}

		for (int j = 0; j < 3; ++j)
		{
			tmp.erase(ii[j]);
			m_faceRingIndices[i].push_back(ii[j]);
		}

		for (auto iter = tmp.begin(); iter != tmp.end(); ++iter)
		{
			m_faceRingIndices[i].push_back(*iter);
		}
	}


}

void TetMesh::computeANs(int tetid)
{
	int edge[3][4] = {{1,0,2,3},{2,0,3,1},{3,0,1,2}};
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
	}
}

Eigen::Matrix3d TetMesh::computeDeformationGradient(int i) // i is tet's index
{
	Eigen::Matrix3d Ds, tmp;
	Ds.col(0) = m_nodes.col(m_tets(1, i)) - m_nodes.col(m_tets(0, i));
	Ds.col(1) = m_nodes.col(m_tets(2, i)) - m_nodes.col(m_tets(0, i));
	Ds.col(2) = m_nodes.col(m_tets(3, i)) - m_nodes.col(m_tets(0, i));

	tmp = Ds * m_Dm_inverses.block<3, 3>(0, 3 * i);

	//for (int k = 0; k < 3; ++k)
	//{
	//	for (int j = 0; j < 3; ++j)
	//	{
	//		tmp(k, j) = fixPrecision(tmp(k, j));
	//	}
	//}

	return tmp;
}


void TetMesh::computeForces()
{
	m_nodes_forces = m_nodes_gravity;
	Eigen::Matrix3d Ds;
	Eigen::Matrix3d F;									// Ds * Dm_inverse
	Eigen::Matrix3d G;									// Green Strain 1/2 * (F_tanspose * F - I)
	Eigen::Matrix3d sigma = Eigen::Matrix3d::Zero();	// Cauch stress using isotropic linear elastic material
	Eigen::Matrix3d P;									// First Piola-Kirchhoff stress
	for (int i = 0; i < n_tets; ++i)
	{
		Ds.col(0) = m_nodes.col(m_tets(1, i)) - m_nodes.col(m_tets(0, i));
		Ds.col(1) = m_nodes.col(m_tets(2, i)) - m_nodes.col(m_tets(0, i));
		Ds.col(2) = m_nodes.col(m_tets(3, i)) - m_nodes.col(m_tets(0, i));
		
		F = Ds * m_Dm_inverses.block<3,3>(0,3*i);
		G = (F.transpose()*F - Eigen::Matrix3d::Identity()) / 2.0;
		//C = (F.transpose() + F) / 2.0 - Eigen::Matrix3d::Identity();
		//test
		//std::cout << "green strain: " << G << std::endl;

		sigma(0, 0) = G(0, 0) *m_Enu3 + G(1, 1)*m_Enu2 + G(2, 2)*m_Enu2;
		sigma(1, 1) = G(0, 0) *m_Enu2 + G(1, 1)*m_Enu3 + G(2, 2)*m_Enu2;
		sigma(2, 2) = G(0, 0) *m_Enu2 + G(1, 1)*m_Enu2 + G(2, 2)*m_Enu3;

		sigma(1, 2) = sigma(2, 1) = G(1, 2) * m_Enu1;
		sigma(0, 2) = sigma(2, 0) = G(0, 2) * m_Enu1;
		sigma(0, 1) = sigma(1, 0) = G(0, 1) * m_Enu1;

		// Cauchy stress to First Piola-Kirchhoff stress
		P = F.determinant() * sigma * F.transpose().inverse();

		//test
		//std::cout << "det(F): " << F.determinant() << std::endl;
		//std::cout << " sigma: " << sigma<<std::endl;
		//std::cout << "F-T: " << F.transpose().inverse() << std::endl;

		// compute forces for node1, node2 and node3.  f4 = - (f1+f2+f3)
		Eigen::Matrix3d inner_forces = P * m_ANs.block<3, 3>(0, 3 * i);

		m_nodes_forces.col(m_tets(1, i)) += inner_forces.col(0);
		m_nodes_forces.col(m_tets(2, i)) += inner_forces.col(1);
		m_nodes_forces.col(m_tets(3, i)) += inner_forces.col(2);
		m_nodes_forces.col(m_tets(0, i)) -= (inner_forces.col(0)+inner_forces.col(1)+inner_forces.col(2));
	}
	//std::cout << "one iteration" << std::endl;
}

void TetMesh::setTetMaterial(double e, double nu, double den)
{
	m_E  = e;
	m_nu = nu;
	m_density = den;

	m_Enu1 = m_E / (1 + m_nu);
	m_Enu2 = m_Enu1 * m_nu / (1- m_nu * 2);
	m_Enu3 = m_Enu1 * (1 - m_nu) / (1 - m_nu * 2);
}

void TetMesh::computeBoundfaceNormals()
{
	Eigen::Vector3d v1, v2, v3;
	for (int i = 0; i < n_bound_faces;++i)
	{
		v1 = m_nodes.col(m_bound_faces(1, i)) - m_nodes.col(m_bound_faces(0, i));
		v2 = m_nodes.col(m_bound_faces(2, i)) - m_nodes.col(m_bound_faces(0, i));
		m_bound_normals.col(i) = v2.cross(v1).normalized();

		m_face_centers.col(i) = (m_nodes.col(m_bound_faces(0, i))
								+m_nodes.col(m_bound_faces(1, i))
								+m_nodes.col(m_bound_faces(2, i)))/3.0;
	}
}

int TetMesh::pickFacebyRay(const Eigen::Vector3d &orig, const Eigen::Vector3d &direct)
{
	double max_distance = 1e8;
	double r = (m_nodes.col(m_bound_faces(1, 0)) - m_nodes.col(m_bound_faces(0, 0))).norm();
	double r2 = r*r;
	double d2,a,b;
	Eigen::Vector3d o_v;
	int pickID = -1;

	for (int i = 0; i < n_bound_faces; ++i)
	{
		if (m_bound_normals.col(i).dot(direct) < 0)
		{
			Eigen::Vector3d tmp = m_face_centers.col(i);
			o_v = orig - m_face_centers.col(i);
			a = o_v.dot(o_v);
			b = o_v.dot(direct);
			d2 = a - b*b;
			if (d2 < r2 && d2 < max_distance)
			{
				pickID = i;
				max_distance = d2;
			}
		}
	}
	return pickID;
}

void TetMesh::dragFace(int faceID, const Eigen::Vector3d &dragline)
{
	double r = (m_nodes.col(m_bound_faces(1, 0)) - m_nodes.col(m_bound_faces(0, 0))).norm();
	double factor = dragline.norm()/r * 100;
	if (factor > 1000)
		factor = 1000;
	//std::cout << factor << std::endl;
	Eigen::Vector3d dir = dragline.normalized();
	//std::cout << dir << std::endl;
	//if (factor > 10)
	//	factor = 10;

	resetExternalForce();
	int nodeID;
	for (int i = 0; i < 3; ++i)
	{
		nodeID = m_bound_faces.col(faceID)[i];
		//std::cout << nodeID << std::endl;
		m_nodes_external_forces.col(nodeID) += factor * m_nodes_mass(nodeID) * dir;
	}
}

void TetMesh::dragFaceRing(int faceID, const Eigen::Vector3d &dragline)
{
	double r = (m_nodes.col(m_bound_faces(1, 0)) - m_nodes.col(m_bound_faces(0, 0))).norm();
	double factor = dragline.norm() / r * 100;
	if (factor > 800)
		factor = 800;
	//std::cout << factor << std::endl;
	Eigen::Vector3d dir = dragline.normalized();
	//std::cout << dir << std::endl;
	//if (factor > 10)
	//	factor = 10;

	resetExternalForce();
	int nodeID;
	for (int i = 0; i < 3; ++i)
	{
		nodeID = m_faceRingIndices[faceID][i];
		//std::cout << nodeID << std::endl;
		m_nodes_external_forces.col(nodeID) += factor * m_nodes_mass(nodeID) * dir;
	}

	for (auto iter = m_faceRingIndices[faceID].begin() + 3; iter != m_faceRingIndices[faceID].end(); ++iter)
	{
		m_nodes_external_forces.col(*iter) += 0.5 * factor * m_nodes_mass(*iter) * dir;
	}
}

void TetMesh::updateNodesVelocities(const Eigen::MatrixXd & pos, const Eigen::MatrixXd & vel)
{
	m_nodes = pos;
	m_velocities = vel;

	//cube test with fixed node 
	//m_nodes.col(8) = Eigen::Vector3d(0, 0.5, -0.5);
	//m_velocities.col(8) = Eigen::Vector3d::Zero();

	// tet test with fixed node
	//m_nodes.col(0) = Eigen::Vector3d(1.0, 1.0, 1.0);
	//m_velocities.col(0) = Eigen::Vector3d::Zero();
	//m_nodes.col(1) = Eigen::Vector3d(-1.0, 1.0, -1.0);
	//m_velocities.col(1) = Eigen::Vector3d::Zero();
}

void TetMesh::drawTetBoundFace()
{
	computeBoundfaceNormals();

	glColor3f(0.251, 0.424, 0.7);
	for (int i = 0; i < n_bound_faces; ++i)
	{
		
		glBegin(GL_TRIANGLES);
		glNormal3dv(m_bound_normals.col(i).data());
		glVertex3dv(m_nodes.col(m_bound_faces(0, i)).data());
		glVertex3dv(m_nodes.col(m_bound_faces(1, i)).data());
		glVertex3dv(m_nodes.col(m_bound_faces(2, i)).data());
		glEnd();
	}
}

void TetMesh::drawDraggedNodes(int faceID)
{
	std::vector<int> &ids = m_faceRingIndices[faceID];
	glColor3f(1.00,0.25,1.00);
	glPointSize(2.0);
	for (auto iter = ids.begin(); iter != ids.end(); ++iter)
	{
		glBegin(GL_POINTS);
		glVertex3dv(m_nodes.col(*iter).data());
		glEnd();
	}
}

void TetMesh::writeMatrix(QString filename, Eigen::MatrixXd mat)
{
	QFile file(filename);
	if (file.open(QIODevice::ReadWrite)) {
		QTextStream stream(&file);
		int row = mat.rows();
		int col = mat.cols();

		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
			{
				stream << mat(i, j) << ",";
			}
			stream << endl;
		}
		file.close();
	}
}

double TetMesh::fixPrecision(double m)
{
	bool isNegtive = false;
	if (m < 0)
	{
		isNegtive = true;
		m *= -1.0;
	}

	double itg = std::floor(m);
	double f = m - itg;
	f = std::round(f*1e8) * 1e-8;
	itg += f;

	if (isNegtive)
		return -itg;
	else
		return itg;
}

void TetMesh::writeNodes(QString nodefile)
{
	QFile file(nodefile);
	if (file.open(QIODevice::ReadWrite)) {
		QTextStream stream(&file);
		stream << n_nodes << " 3 0 0" << endl;
		for (size_t i = 0; i < n_nodes; ++i)
		{
			stream << i << " " << m_nodes(0, i) << " " << m_nodes(1, i) << " " << m_nodes(2, i) << endl;
		}
		file.close();
	}
}

