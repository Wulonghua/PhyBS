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

	Eigen::VectorXf projectLocalConstraints(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
		const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
		const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext, int num_threads);

	//collect each contraint's projection to a large memory and then sum up with each node. This way avoids atomic operation in parallelism
	Eigen::VectorXf projectLocalConstraints2(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
		const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
		const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext);

	Eigen::VectorXf projectLocalConstraints2(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
		const Eigen::MatrixXi &tets, float t, Eigen::MatrixXf s, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
		const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext, int num_threads);

	void solveGlobalStep(Eigen::MatrixXf &pos, Eigen::VectorXf &b);


	void initLocalProjection(const Eigen::MatrixXi & tets);

	Eigen::SparseMatrix<float, Eigen::RowMajor> m_globalSolverMat;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<float, Eigen::RowMajor>> m_pardiso_solver;
	Eigen::VectorXf m_stiffWeight;
	Eigen::MatrixXf m_pos_new;
	Eigen::MatrixXf m_localProjections;  // store local projections result for each tet
	std::vector<std::vector<int>> m_proj_idx; // each node's projection index in the plain array of m_localProjections. 

	float m_stiffness;

	int n_nodes;
	int n_tets;
	int m_iterations;

	//for test
	QElapsedTimer m_time;
	int m_elapses1;
	int m_elapses2;
	
};

