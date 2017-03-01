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

	void reset();
	void doDescentOpt(int iterations);

private:
	Eigen::MatrixXf computeGradient(const Eigen::MatrixXf &ext_f);
	void initialization();
	float computeTotalEnergy(const Eigen::MatrixXf &pos);

	Eigen::MatrixXf m_posk_next;
	Eigen::MatrixXf m_posk;
	Eigen::MatrixXf m_posk_last;
	Eigen::MatrixXf m_pos;
	Eigen::MatrixXf m_posk_old;
	//Eigen::MatrixXf m_pos0;
	//Eigen::MatrixXf m_pos1;

	Eigen::MatrixXf m_vel_old;

	float m_h;

	std::shared_ptr<TetMesh>					 m_tetMesh;
	std::shared_ptr<IsotropicNeohookeanMaterial> m_IsoMaterial;

	int n_nodes;
	int n_tets;

	float m_alpha;

	float m_energy;
	float m_energy_old;

	float m_rho;
	float m_omega;

	int		m_profile_k[3];
	float	m_profile_v[3];
	int		m_profile_n;
};

