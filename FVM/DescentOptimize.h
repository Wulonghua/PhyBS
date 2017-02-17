#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/svd>
#include <vector>
#include <memory>

#include "TetMesh.h"
#include "IsotropicNeohookeanMaterial.h"

class DescentOptimize
{
public:
	DescentOptimize(std::shared_ptr<TetMesh> tetMesh, std::shared_ptr<IsotropicNeohookeanMaterial> isoMaterial);
	~DescentOptimize();

	void initialization();


private:
	Eigen::MatrixXf computeGradient();

	Eigen::MatrixXf m_posk;
	//Eigen::MatrixXf m_pos0;
	//Eigen::MatrixXf m_pos1;

	Eigen::MatrixXf m_vel_old;

	float m_timeStep;

	std::shared_ptr<TetMesh>					 m_tetMesh;
	std::shared_ptr<IsotropicNeohookeanMaterial> m_IsoMaterial;

};

