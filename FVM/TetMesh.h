#pragma once

#include "Eigen/Dense"
#include "qfile.h"
#include "qstring.h"
#include "qtextstream.h"
#include "qopengl.h"

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

private:
	void initModel();
	void computeANs(int tetid);
	void computeBoundfaceNormals();

	Eigen::MatrixXd m_nodes;			// nodes' positions     : 3*n matrix
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
