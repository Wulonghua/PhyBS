#pragma once
#include "Eigen/Dense"

class TimeIntegration
{
public:
	TimeIntegration(int n_nodes);
	~TimeIntegration();

	void simuExplicit(	const Eigen::MatrixXd & pos,
						const Eigen::MatrixXd & vel,
						const Eigen::MatrixXd & force,
						const Eigen::VectorXd & mass);
	void setTimeStep(double t) { m_t = t; }
	Eigen::MatrixXd & getPositions() { return m_positions;}
	Eigen::MatrixXd & getVelocities() { return m_velocities; }
private:
	Eigen::MatrixXd m_positions;  // nodes' positions after one time step
	Eigen::MatrixXd m_velocities; // nodes' velocities after one time step  

	double m_t;                   // time step
	int n_nodes;
};