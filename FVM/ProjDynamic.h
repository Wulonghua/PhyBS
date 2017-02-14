#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/svd>
#include <Eigen/PardisoSupport>
#include <vector>
#include <memory>
#include <qelapsedtimer.h>

#include "TetMesh.h"

class ProjDynamic
{
public:
	ProjDynamic(std::shared_ptr<TetMesh> tetMesh);
	virtual ~ProjDynamic();

	void buildGlobalSolverMatrix(const Eigen::VectorXf & node_mass, const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Bm);
	void doProjDynamics(Eigen::MatrixXf &pos, Eigen::MatrixXf &vel,
						const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
						const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Dm_inverse, const Eigen::MatrixXf & fext);
private:

	Eigen::VectorXf projectLocalConstraints(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
		const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
		const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext);
	void solveGlobalStep(Eigen::MatrixXf &pos, Eigen::VectorXf &b);

	Eigen::SparseMatrix<float, Eigen::RowMajor> m_globalSolverMat;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<float, Eigen::RowMajor>> m_pardiso_solver;
	Eigen::VectorXf m_stiffWeight;
	Eigen::MatrixXf m_pos_new;
	float m_stiffness;

	int n_nodes;
	int n_tets;
	int m_iterations;


	//for test
	QElapsedTimer m_time;
	int m_elapses1;
	int m_elapses2;
	
};

