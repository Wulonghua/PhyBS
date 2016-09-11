#include "TimeIntegration.h"


TimeIntegration::TimeIntegration(int num_nodes) : m_t(0.001)
{
	m_positions = Eigen::MatrixXd::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXd::Zero(3, num_nodes);
	n_nodes = num_nodes;
}


TimeIntegration::~TimeIntegration()
{
}

void TimeIntegration::simuExplicit(	const Eigen::MatrixXd & pos,
									const Eigen::MatrixXd & vel,
									const Eigen::MatrixXd & force,
									const Eigen::VectorXd & mass)
{
	for (int i = 0; i < n_nodes; ++i)
	{
		m_velocities.col(i) = vel.col(i) + m_t * force.col(i) / mass(i);
		m_positions.col(i) = pos.col(i) + m_velocities.col(i) * m_t;
	}
}