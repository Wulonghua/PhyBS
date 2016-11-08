#pragma once

#include <Eigen/Dense>
#include <qfile.h>
#include <qstring.h>
#include <qtextstream.h>
#include <qopengl.h>
#include <qfile.h>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

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

	void setTetMaterial(float e, float nu, float den);

	// using the simplest linear isotropic model
	void computeForces();
	Eigen::Matrix3f computeDeformationGradient(int tetID);

	void drawTetBoundFace();
	void drawDraggedNodes(int faceID);
	void updateNodesVelocities(const Eigen::MatrixXf & pos, const Eigen::MatrixXf & vel);

	// The following are inline functions
	Eigen::MatrixXf & getNodes() { return m_nodes; }
	Eigen::MatrixXf & getVelocities() { return m_velocities; }
	Eigen::MatrixXf & getForces() { return m_nodes_forces; }
	Eigen::VectorXf & getMasses() { return m_nodes_mass; }
	Eigen::MatrixXf & getRestPosition() { return m_rest_positions; }
	std::vector<int> & getConstraintIDs() { return m_constraintIDs; }
	int getNodesNum() { return n_nodes; }
	int getTetsNum() { return n_tets; }
	float getE() { return m_E; }
	float getNu() { return m_nu; }
	Eigen::Vector3f getFaceCenter(int i) { return m_face_centers.col(i); }
	Eigen::Matrix3f getAN(int tetID) { return m_ANs.block<3, 3>(0, 3 * tetID); }
	Eigen::Matrix3f getDmInv(int tetID) { return m_Dm_inverses.block<3, 3>(0, 3 * tetID); }
	int getNodeGlobalIDinTet(int tetID, int localID) { return m_tets(localID, tetID); }
	

	void addNodesForces(Eigen::MatrixXf &forces) { m_nodes_forces += forces;}
	void addNodeForce(int nodeID, const Eigen::Vector3f &force) { m_nodes_forces.col(nodeID) += force; }
	void initForcesFromGravityExternals() { m_nodes_forces = m_nodes_gravity + m_nodes_external_forces; }
	void resetExternalForce() { m_nodes_external_forces.setZero();}
	void dragFace(int faceID, const Eigen::Vector3f &dragline);
	void dragFaceRing(int faceID, const Eigen::Vector3f &dragline);
	// approximate method, shoot the ray and test the center of each face, if it is within certain distance
	// then the face is chosen.
	int pickFacebyRay(const Eigen::Vector3f &orig, const Eigen::Vector3f &direct);

	// for test
	void writeMatrix(QString file, Eigen::MatrixXf mat);
	void writeNodes(QString file);

	float fixPrecision(float m);

private:
	void initModel();
	void computeANs(int tetid);
	void computeBoundfaceNormals();
	void computeBoundfaceRingIndices();
	

	Eigen::MatrixXf m_nodes;			// nodes' positions     : 3*n matrix
	Eigen::MatrixXf m_rest_positions;	// nodes' restposition	: 3*n matrix
	Eigen::MatrixXf m_velocities;		// nodes' velocities    : 3*n matrix
	Eigen::MatrixXi m_tets;				// tetrahedra's indices : 4*m matrix

	Eigen::MatrixXi m_bound_faces;		// the surfaces' indices of tet-mesh's boundary.
	Eigen::MatrixXf m_bound_normals;    // the surfaces' normals.
	Eigen::MatrixXf m_face_centers;
	std::vector<std::vector<int>> m_faceRingIndices;   // store the indices of the nodes which connect to certain bound face.

	Eigen::MatrixXf m_Dm_inverses;
	Eigen::MatrixXf m_ANs;

	Eigen::VectorXf m_nodes_mass;		// nodes 
	Eigen::MatrixXf m_nodes_gravity;
	Eigen::MatrixXf m_nodes_external_forces;

	Eigen::MatrixXf m_nodes_forces;
	
	std::vector<int> m_constraintIDs;

	float m_density;
	float m_E;
	float m_nu;

	float m_Enu1, m_Enu2, m_Enu3;    // for linear 

	int n_nodes;
	int n_tets;
	int n_bound_faces;

};
