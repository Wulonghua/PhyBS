#include "TimeIntegration.h"


TimeIntegration::TimeIntegration(int num_nodes) : m_t(1e-6)
{
	m_positions = Eigen::MatrixXd::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXd::Zero(3, num_nodes);
	n_nodes = num_nodes;
}

TimeIntegration::TimeIntegration(int num_nodes, Eigen::VectorXd m) : m_t(1e-3)
{
	m_positions = Eigen::MatrixXd::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXd::Zero(3, num_nodes);
	n_nodes = num_nodes;
	m_masses = Eigen::VectorXd::Zero(3 * num_nodes);

	for (size_t i = 0; i < num_nodes; ++i)
	{
		m_masses[3 * i] = m[i];
		m_masses[3 * i + 1] = m[i];
		m_masses[3 * i + 2] = m[i];
	}

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

void TimeIntegration::BackEuler( Eigen::MatrixXd & pos,
								 Eigen::MatrixXd & vel,
								 Eigen::MatrixXd & force,
								 Eigen::MatrixXd & K)
{
	int n = pos.cols();
	Eigen::MatrixXd A = -m_t*m_t*K;
	A.diagonal() += m_masses;
	//std::cout << A << std::endl;
	pos.resize(3*n,1);
	vel.resize(3*n,1);
	force.resize(3*n,1);
	Eigen::VectorXd p = pos;
	Eigen::VectorXd v = vel;
	Eigen::VectorXd f = force;

	Eigen::VectorXd b = m_masses.cwiseProduct(v) + m_t * f;

	Eigen::ConjugateGradient<Eigen::MatrixXd> cg_solver;
	cg_solver.compute(A);
	cg_solver.setTolerance(1e-8);
	v = cg_solver.solve(b);
	p += m_t * v;

	pos = p.matrix();
	vel = v.matrix();
	pos.resize(3, n);
	vel.resize(3 ,n);
	force.resize(3,n);
}