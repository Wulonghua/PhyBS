#pragma once

#include "Eigen/Dense"
#include <iostream>
#include <qfile.h>
#include <qstring.h>
#include <qtextstream.h>

class TetMesh
{
public:
	TetMesh();
	~TetMesh();
	void LoadNodesFromFile(QString filename);
	void LoadTetsFromFile(QString filename);
	void LoadFacesFromFile(QString filename);
	void SetTetMaterial(double e, double nu);
	void ComputeForces();



private:
	void InitModel();
	void ComputeANs(int tetid);

	Eigen::MatrixXd m_nodes;			// nodes' positions     : 3*n matrix
	Eigen::MatrixXi m_tets;				// tetrahedra's indices : 4*m matrix
	Eigen::MatrixXi m_bound_faces;		// the surfaces' indices of tet-mesh's boundary.

	Eigen::MatrixXd m_Dm_inverses;
	Eigen::MatrixXd m_ANs;

	Eigen::VectorXd m_nodes_mass;
	Eigen::MatrixXd m_nodes_gravity;

	Eigen::MatrixXd m_nodes_forces;


	double m_E;
	double m_nu;

	int n_nodes;
	int n_tets;
	int n_bound_faces;
};

