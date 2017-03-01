#include "DescentOptimize.h"


DescentOptimize::DescentOptimize(std::shared_ptr<TetMesh> tetMesh, std::shared_ptr<IsotropicNeohookeanMaterial> isoMaterial):
m_h(0.03), m_energy(0.0), m_energy_old(0.0)
{
	m_tetMesh = tetMesh;
	m_IsoMaterial = isoMaterial;

	n_nodes = m_tetMesh->getNodesNum();
	n_tets = m_tetMesh->getTetsNum();

	m_profile_k[0] = 1;
	m_profile_k[1] = 7;
	m_profile_k[2] = 10;
	m_profile_v[0] = 0.96;
	m_profile_v[1] = 0.991;
	m_profile_v[2] = 0.999;
	m_profile_n = 3;

	reset();
}


DescentOptimize::~DescentOptimize()
{
}

void DescentOptimize::reset()
{
	m_vel_old = m_tetMesh->getVelocities();
	m_posk_last = m_posk = m_tetMesh->getNodes();
}

void DescentOptimize::initialization()
{
	//m_posk = m_pos1 + m_timeStep * m_vel1 + 0.2 * m_timeStep * (m_vel1 - m_vel0);
	m_pos = m_tetMesh->getNodes() + m_h * m_tetMesh->getVelocities();
	m_posk = m_pos + 0.2 * m_h* (m_tetMesh->getVelocities()-m_vel_old);
	//m_posk_old = m_posk;
}

Eigen::MatrixXf DescentOptimize::computeGradient(const Eigen::MatrixXf &ext_f)
{
	Eigen::MatrixXf g = m_posk - m_pos;
	for (int i = 0; i < n_nodes; ++i)
	{
		g.col(i) = g.col(i) * m_tetMesh->getMasses()[i] / (m_h*m_h);
	}

	Eigen::MatrixXf f = m_IsoMaterial->computeInnerForceFromPos(m_posk);
	g = g - f - ext_f;

	auto consIDs = m_tetMesh->getConstraintIDs();
	for (int i = 0; i < consIDs.size(); ++i)
	{
		g.col(consIDs[i]).setZero();
	}

	return g;
}

float DescentOptimize::computeTotalEnergy(const Eigen::MatrixXf &pos)
{
	float E = 0.0;
	m_IsoMaterial->computeElasticEnergyFromPos(pos);
	Eigen::VectorXf Elastic = m_IsoMaterial->getElasticEnergys();

	// elastic energy
	for (int i = 0; i < n_tets; ++i)
		E += Elastic[i];

	// kinetic energy +  gravitational potential energy
	Eigen::MatrixXf p = pos - m_pos;
	auto consIDs = m_tetMesh->getConstraintIDs();
	int fixi = 0;
	float h2_inv = 1 / (m_h*m_h);
	for (int i = 0; i < n_nodes; ++i)
	{
		if (consIDs.size() > 0 && consIDs[fixi] == i)
		{
			fixi++;
			continue;
		}
		E += 0.5 * m_tetMesh->getMasses()[i] * ( p.col(i).squaredNorm() * h2_inv + 2 * 9.8 * pos(1,i));
	}
	E += m_tetMesh->computeExternalForcesEnergy(pos);
	return E;
}

void DescentOptimize::doDescentOpt(int iterations)
{
	Eigen::MatrixXf delta_q;
	Eigen::VectorXf H_q;
	float m_h2_inv;
	m_alpha = 0.1;
	//m_IsoMaterial->computeElasticEnergyFromPos(m_posk);
	//Eigen::VectorXf tmp = m_IsoMaterial->getElasticEnergys();
	initialization();
	//std::cout << std::endl;
	m_tetMesh->initForcesFromGravityExternals();
	Eigen::MatrixXf extForce = m_tetMesh->getForces();
	for (int k = 0; k < iterations; ++k)
	{
		if (m_alpha < 0.01) break;
		
		if (k % 8 == 0)
		{
			m_energy = computeTotalEnergy(m_posk);
			
			if (k != 0 && m_energy > m_energy_old)
			{
				m_alpha *= 0.7;
				iterations -= k-8;
				k -= 1;
				m_posk = m_posk_old;
				//std::cout << "restore to last state: " << computeTotalEnergy(m_posk) << std::endl;
				//std::cout << "k = " << k << std::endl;
				//std::cout << "alpha =" << m_alpha << std::endl;
				//std::cout << "iterations = " << iterations << std::endl;
				continue;
			}
			else
			{
				m_energy_old = m_energy;
				m_posk_old = m_posk;
			}
		}

		//std::cout << "energy: " << computeTotalEnergy(m_posk) << std::endl;
		if (k % 32 == 0)
		{
			H_q = m_IsoMaterial->computeGlobalStiffnessMatrixFromPos(m_posk).diagonal();
			for (int i = 0; i < n_nodes; ++i)
			{
				m_h2_inv = m_tetMesh->getMasses()[i] / (m_h * m_h);
				H_q[3 * i + 0] += m_h2_inv;
				H_q[3 * i + 1] += m_h2_inv;
				H_q[3 * i + 2] += m_h2_inv;
			}
		}
		delta_q = computeGradient(extForce);
		for (int i = 0; i < n_nodes; ++i)
		{
			delta_q(0, i) /= -H_q[3 * i + 0];
			delta_q(1, i) /= -H_q[3 * i + 1];
			delta_q(2, i) /= -H_q[3 * i + 2];
		}

		m_posk_next = m_posk + m_alpha * delta_q;

		//if (k < 10) m_omega = 1.0;
		//else if(k == 10) m_omega = 2 / (2 - m_rho*m_rho);
		//else m_omega = 4 / (4 - m_rho*m_rho*m_omega);

		if (k == 0)				
		{ 
			m_rho = 0;
			m_omega = 1; 
		}
		m_omega = 4 / (4 - m_rho*m_rho*m_omega);
		for (int i = 0; i<m_profile_n; i++)
		{
			if (k == m_profile_k[i] - 1)
			{
				m_rho = 0;
				m_omega = 1;
			}
			if (k == m_profile_k[i])
			{
				m_rho = m_profile_v[i];
				m_omega = 2 / (2 - m_rho*m_rho);
				break;
			}
		}

		if (m_omega != 1.0)
		m_posk_next = m_omega * (m_posk_next - m_posk_last) + m_posk_last;

		m_posk_last = m_posk;
		m_posk = m_posk_next;
	}
	m_vel_old = m_tetMesh->getVelocities();
	m_tetMesh->getVelocities() = 0.99 * (m_posk - m_tetMesh->getNodes()) / m_h;
	m_tetMesh->getNodes() = m_posk;
	m_tetMesh->resetExternalForce();
	//std::cout << m_posk(1, 0) << std::endl;
}
