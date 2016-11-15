#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/PardisoSupport>
#include <iostream>
#include "mkl.h"

class TimeIntegration
{
public:
	TimeIntegration(int n_nodes);
	TimeIntegration(int n_nodes, Eigen::VectorXf m_masses);
	TimeIntegration(int n_nodes, Eigen::VectorXf m_masses, std::vector<int> constraints, Eigen::MatrixXf rest_pos);
	~TimeIntegration();

	void simuExplicit(const Eigen::MatrixXf & pos,
		const Eigen::MatrixXf & vel,
		const Eigen::MatrixXf & force,
		const Eigen::VectorXf & mass);
	void BackEuler(Eigen::MatrixXf & pos,
		Eigen::MatrixXf & vel,
		Eigen::MatrixXf & force,
		Eigen::MatrixXf & K);
	void BackEuler(Eigen::MatrixXf & pos,
		Eigen::MatrixXf & vel,
		Eigen::MatrixXf & force,
		Eigen::SparseMatrix<float, Eigen::RowMajor> & K);

	void BackEuler(Eigen::MatrixXf & pos,
		Eigen::MatrixXf & restPos,
		Eigen::MatrixXf & vel,
		Eigen::MatrixXf & force,
		Eigen::SparseMatrix<float, Eigen::RowMajor> & K);

	void setTimeStep(float t) { m_t = t; }
	Eigen::MatrixXf  getPositions() { return m_positions; }
	Eigen::MatrixXf  getVelocities() { return m_velocities; }
	float & getTimeStep() { return m_t; }

private:
	void addGroundConstraints(float y, Eigen::MatrixXf & pos, Eigen::MatrixXf & vel);

	Eigen::PardisoLDLT<Eigen::SparseMatrix<float, Eigen::RowMajor>> m_pardiso_solver;

	Eigen::MatrixXf m_positions;  // nodes' positions after one time step
	Eigen::MatrixXf m_velocities; // nodes' velocities after one time step
	Eigen::VectorXf m_masses;

	float m_dumpingAlpha;
	float m_dumpingBelta;

	float m_t;                   // time step
	int n_nodes;

	std::vector<int> m_constraints;
	Eigen::MatrixXf	 m_rest;

};