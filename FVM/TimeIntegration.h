#pragma once

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/PardisoSupport>
#include <iostream>

class TimeIntegration
{
public:
	TimeIntegration(int n_nodes);
	TimeIntegration(int n_nodes, Eigen::VectorXd m_masses);
	TimeIntegration(int n_nodes, Eigen::VectorXd m_masses, std::vector<int> constraints, Eigen::MatrixXd rest_pos);
	~TimeIntegration();

	void simuExplicit(const Eigen::MatrixXd & pos,
		const Eigen::MatrixXd & vel,
		const Eigen::MatrixXd & force,
		const Eigen::VectorXd & mass);
	void BackEuler(Eigen::MatrixXd & pos,
		Eigen::MatrixXd & vel,
		Eigen::MatrixXd & force,
		Eigen::MatrixXd & K);
	void BackEuler(Eigen::MatrixXd & pos,
		Eigen::MatrixXd & vel,
		Eigen::MatrixXd & force,
		Eigen::SparseMatrix<double> & K);

	void BackEuler(Eigen::MatrixXd & pos,
		Eigen::MatrixXd & restPos,
		Eigen::MatrixXd & vel,
		Eigen::MatrixXd & force,
		Eigen::SparseMatrix<double> & K);

	void setTimeStep(double t) { m_t = t; }
	Eigen::MatrixXd  getPositions() { return m_positions; }
	Eigen::MatrixXd  getVelocities() { return m_velocities; }
	double & getTimeStep() { return m_t; }

private:
	void addGroundConstraints(double y, Eigen::MatrixXd & pos, Eigen::MatrixXd & vel);

	Eigen::MatrixXd m_positions;  // nodes' positions after one time step
	Eigen::MatrixXd m_velocities; // nodes' velocities after one time step
	Eigen::VectorXd m_masses;

	double m_dumpingAlpha;
	double m_dumpingBelta;

	double m_t;                   // time step
	int n_nodes;

	std::vector<int> m_constraints;
	Eigen::MatrixXd	 m_rest;

};