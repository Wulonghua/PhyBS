#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/svd>
#include <Eigen/PardisoSupport>
#include <vector>
#include <memory>

#include "TetMesh.h"

class ProjDynamic
{
public:
	ProjDynamic(std::shared_ptr<TetMesh> tetMesh);
	virtual ~ProjDynamic();

	void buildGlobalSolverMatrix(const Eigen::VectorXf & node_mass, const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Bm);
	Eigen::VectorXf projectLocalConstraints(const Eigen::VectorXf & node_mass, const Eigen::VectorXf &node_inv_mass,
								const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &pos, const Eigen::MatrixXf &Dm_inverse,
								const Eigen::MatrixXf &vel, const Eigen::MatrixXf & fext);
	void solveGlobalStep(Eigen::MatrixXf &pos, const Eigen::VectorXf &b);
private:

	Eigen::SparseMatrix<float, Eigen::RowMajor> m_globalSolverMat;
	Eigen::PardisoLDLT<Eigen::SparseMatrix<float, Eigen::RowMajor>> m_pardiso_solver;
	Eigen::VectorXf m_stiffWeight;
	float m_stiffness;

	int n_nodes;
	int n_tets;
	
};

