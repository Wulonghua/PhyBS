#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <memory>

#include "TetMesh.h"

class ProjDynamic
{
public:
	ProjDynamic(std::shared_ptr<TetMesh> tetMesh);
	virtual ~ProjDynamic();

	void buildGlobalSolverMatrix(const Eigen::VectorXf & node_mass, const Eigen::MatrixXi &tets, float t, const Eigen::MatrixXf &Bm);
private:

	Eigen::SparseMatrix<float, Eigen::RowMajor> m_globalSolverMat;
	Eigen::VectorXf m_stiffWeight;
	float m_stiffness;

	int n_nodes;
	int n_tets;
	
};

