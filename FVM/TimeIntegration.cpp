#include "TimeIntegration.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

TimeIntegration::TimeIntegration(int num_nodes) : m_t(1e-6)
{
	m_positions = Eigen::MatrixXf::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXf::Zero(3, num_nodes);
	n_nodes = num_nodes;
}

TimeIntegration::TimeIntegration(int num_nodes, Eigen::VectorXf m) : m_t(2e-3)
{
	m_positions = Eigen::MatrixXf::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXf::Zero(3, num_nodes);
	n_nodes = num_nodes;
	m_masses = Eigen::VectorXf::Zero(3 * num_nodes);

	for (size_t i = 0; i < num_nodes; ++i)
	{
		m_masses[3 * i] = m[i];
		m_masses[3 * i + 1] = m[i];
		m_masses[3 * i + 2] = m[i];
	}

}

TimeIntegration::TimeIntegration(int num_nodes, Eigen::VectorXf m, std::vector<int> constraints, Eigen::MatrixXf rest_pos) :
m_t(0.03), m_constraints(constraints), m_rest(rest_pos), m_dumpingAlpha(0.02), m_dumpingBelta(0.02)
{
	m_positions = Eigen::MatrixXf::Zero(3, num_nodes);
	m_velocities = Eigen::MatrixXf::Zero(3, num_nodes);
	n_nodes = num_nodes;
	m_masses = Eigen::VectorXf::Zero(3 * num_nodes);

	for (size_t i = 0; i < num_nodes; ++i)
	{
		m_masses[3 * i] = m[i];
		m_masses[3 * i + 1] = m[i];
		m_masses[3 * i + 2] = m[i];
	}

	for (size_t i = 0; i < m_constraints.size(); ++i)
	{
		m_masses[3 * m_constraints[i]] = 1e8;
		m_masses[3 * m_constraints[i] + 1] = 1e8;
		m_masses[3 * m_constraints[i] + 2] = 1e8;
	}
}


TimeIntegration::~TimeIntegration()
{
}

void TimeIntegration::simuExplicit(	const Eigen::MatrixXf & pos,
									const Eigen::MatrixXf & vel,
									const Eigen::MatrixXf & force,
									const Eigen::VectorXf & mass)
{
	for (int i = 0; i < n_nodes; ++i)
	{
		m_velocities.col(i) = vel.col(i) + m_t * force.col(i) / mass(i);
		m_positions.col(i) = pos.col(i) + m_velocities.col(i) * m_t;
	}
}

void TimeIntegration::BackEuler(Eigen::MatrixXf & pos,
	Eigen::MatrixXf & vel,
	Eigen::MatrixXf & force,
	Eigen::SparseMatrix<float, Eigen::RowMajor> & K)
{

	int n = pos.cols();

	//std::cout << "K:" << std::endl;
	//std::cout << K << std::endl;

	Eigen::SparseMatrix<float> A = (m_dumpingBelta + m_t)*m_t*K;
	for (int i = 0; i < 3 * n; ++i)
	{
		A.coeffRef(i, i) += (1 + m_t*m_dumpingAlpha)*m_masses(i);
	}

	//std::cout << "LHS:" << std::endl;
	//std::cout << A << std::endl;

	Eigen::Map<Eigen::VectorXf> p(pos.data(), n_nodes * 3);
	Eigen::Map<Eigen::VectorXf> v(vel.data(), n_nodes * 3);
	Eigen::Map<Eigen::VectorXf> f(force.data(), n_nodes * 3);

	//std::cout << "forces" << std::endl;
	//std::cout << f << std::endl;

	Eigen::VectorXf b = m_masses.cwiseProduct(v) + m_t * f;
	
	//for (int i = 0; i < 3; ++i)
	//	printf("%.10f %.10f\n", v[i], b[i]);

	//std::cout << "RHS:" << std::endl;
	//std::cout << b << std::endl;
	
	mkl_set_dynamic(0);
	mkl_set_num_threads(4);
	m_pardiso_solver.compute(A);
	v = m_pardiso_solver.solve(b);

	//std::cout << "v:" << std::endl;
	//std::cout << v << std::endl;

	p += m_t * v;

	//std::cout << "p:" << std::endl;
	//std::cout << p << std::endl;
	
	//std::cout << mkl_get_max_threads() << std::endl;
}

void TimeIntegration::addGroundConstraints(float y, Eigen::MatrixXf & pos, Eigen::MatrixXf & vel)
{
	for (size_t i = 0; i < n_nodes; ++i)
	{
		if (pos(1, i) < y)
		{
			pos(1, i) = y;
			vel.col(i) = Eigen::Vector3f::Zero();
		}
	}
}