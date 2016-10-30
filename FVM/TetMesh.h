#pragma once

#include <Eigen/Dense>
#include <qfile.h>
#include <qstring.h>
#include <qtextstream.h>
#include <qopengl.h>
#include <qfile.h>
#include <iostream>

class RenderWidget;
class TetMesh
{
public:
	TetMesh();
	~TetMesh();
	void initNodesFromFile(QString filename);
	void initTetsFromFile(QString filename);
	void initFacesFromFile(QString filename);
	void updateNodesFromFile(QString filename);

	void setTetMaterial(double e, double nu, double den);

	// using the simplest linear isotropic model
	void computeForces();
	Eigen::Matrix3d computeDeformationGradient(int tetID);

	void drawTetBoundFace();
	void updateNodesVelocities(const Eigen::MatrixXd & pos, const Eigen::MatrixXd & vel);

	// The following are inline functions
	Eigen::MatrixXd & getNodes() { return m_nodes; }
	Eigen::MatrixXd & getVelocities() { return m_velocities; }
	Eigen::MatrixXd & getForces() { return m_nodes_forces; }
	Eigen::VectorXd & getMasses() { return m_nodes_mass; }
	Eigen::MatrixXd & getRestPosition() { return m_rest_positions; }
	std::vector<int> & getConstraintIDs() { return m_constraintIDs; }
	int getNodesNum() { return n_nodes; }
	int getTetsNum() { return n_tets; }
	double getE() { return m_E; }
	double getNu() { return m_nu; }
	Eigen::Vector3d getFaceCenter(int i) { return m_face_centers.col(i); }
	Eigen::Matrix3d getAN(int tetID) { return m_ANs.block<3, 3>(0, 3 * tetID); }
	Eigen::Matrix3d getDmInv(int tetID) { return m_Dm_inverses.block<3, 3>(0, 3 * tetID); }
	int getNodeGlobalIDinTet(int tetID, int localID) { return m_tets(localID, tetID); }
	

	void addNodesForces(Eigen::MatrixXd &forces) { m_nodes_forces += forces;}
	void addNodeForce(int nodeID, Eigen::Vector3d const &force) { m_nodes_forces.col(nodeID) += force; }
	void initForcesFromGravity() { m_nodes_forces = m_nodes_gravity; }

	// approximate method, shoot the ray and test the center of each face, if it is within certain distance
	// then the face is chosen.
	int pickFacebyRay(const Eigen::Vector3d &orig, const Eigen::Vector3d &direct);

	// for test
	void writeMatrix(QString file, Eigen::MatrixXd mat);
	void writeNodes(QString file);

	double fixPrecision(double m);

private:
	void initModel();
	void computeANs(int tetid);
	void computeBoundfaceNormals();
	

	Eigen::MatrixXd m_nodes;			// nodes' positions     : 3*n matrix
	Eigen::MatrixXd m_rest_positions;	// nodes' restposition	: 3*n matrix
	Eigen::MatrixXd m_velocities;		// nodes' velocities    : 3*n matrix
	Eigen::MatrixXi m_tets;				// tetrahedra's indices : 4*m matrix

	Eigen::MatrixXi m_bound_faces;		// the surfaces' indices of tet-mesh's boundary.
	Eigen::MatrixXd m_bound_normals;    // the surfaces' normals.
	Eigen::MatrixXd m_face_centers;

	Eigen::MatrixXd m_Dm_inverses;
	Eigen::MatrixXd m_ANs;

	Eigen::VectorXd m_nodes_mass;		// nodes 
	Eigen::MatrixXd m_nodes_gravity;

	Eigen::MatrixXd m_nodes_forces;
	
	std::vector<int> m_constraintIDs;

	double m_density;
	double m_E;
	double m_nu;

	double m_Enu1, m_Enu2, m_Enu3;    // for linear 

	int n_nodes;
	int n_tets;
	int n_bound_faces;

};
