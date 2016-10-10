#pragma once
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

class TimeIntegration
{
public:
	TimeIntegration(int n_nodes);
	TimeIntegration(int n_nodes, Eigen::VectorXd m_masses);
	~TimeIntegration();

	void simuExplicit(const Eigen::MatrixXd & pos,
		const Eigen::MatrixXd & vel,
		const Eigen::MatrixXd & force,
		const Eigen::VectorXd & mass);
	void BackEuler(Eigen::MatrixXd & pos,
		Eigen::MatrixXd & vel,
		Eigen::MatrixXd & force,
		Eigen::MatrixXd & K);
	void setTimeStep(double t) { m_t = t; }
	Eigen::MatrixXd  getPositions() { return m_positions; }
	Eigen::MatrixXd  getVelocities() { return m_velocities; }
	double & getTimeStep() { return m_t; }
private:
	Eigen::MatrixXd m_positions;  // nodes' positions after one time step
	Eigen::MatrixXd m_velocities; // nodes' velocities after one time step
	Eigen::VectorXd m_masses;

	double m_t;                   // time step
	int n_nodes;
};