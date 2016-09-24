#pragma once

#include <Eigen/Dense>
#include <qfile.h>
#include <qstring.h>
#include <qtextstream.h>
#include <qopengl.h>
#include <iostream>

class RenderWidget;
class TetMesh
{
public:
	TetMesh();
	~TetMesh();
	void loadNodesFromFile(QString filename);
	void loadTetsFromFile(QString filename);
	void loadFacesFromFile(QString filename);
	void setTetMaterial(double e, double nu);
	void computeForces();
	void drawTetBoundFace();
	void updateNodesVelocities(const Eigen::MatrixXd & pos, const Eigen::MatrixXd & vel);
	Eigen::MatrixXd & getNodes() { return m_nodes; }
	Eigen::MatrixXd & getVelocities() { return m_velocities; }
	Eigen::MatrixXd & getForces() { return m_nodes_forces; }
	Eigen::VectorXd & getMasses() { return m_nodes_mass; }
	int getNodesNum() { return n_nodes; }
	int getTetsNum() { return n_tets; }
	double getE() { return m_E; }
	double getNu() { return m_nu; }

private:
	void initModel();
	void computeANs(int tetid);
	void computeBoundfaceNormals();

	Eigen::MatrixXd m_nodes;			// nodes' positions     : 3*n matrix
	Eigen::MatrixXd m_velocities;		// nodes' velocities    : 3*n matrix
	Eigen::MatrixXi m_tets;				// tetrahedra's indices : 4*m matrix

	Eigen::MatrixXi m_bound_faces;		// the surfaces' indices of tet-mesh's boundary.
	Eigen::MatrixXd m_bound_normals;    // the surfaces' normals.


	Eigen::MatrixXd m_Dm_inverses;
	Eigen::MatrixXd m_ANs;

	Eigen::VectorXd m_nodes_mass;
	Eigen::MatrixXd m_nodes_gravity;

	Eigen::MatrixXd m_nodes_forces;


	double m_E;
	double m_nu;

	double m_Enu1, m_Enu2, m_Enu3;    

	int n_nodes;
	int n_tets;
	int n_bound_faces;
};
