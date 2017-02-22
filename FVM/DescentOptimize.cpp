#include "DescentOptimize.h"


DescentOptimize::DescentOptimize(std::shared_ptr<TetMesh> tetMesh, std::shared_ptr<IsotropicNeohookeanMaterial> isoMaterial):
m_iterations(96), m_h(0.03), m_energy(0.0), m_energy_old(0.0), m_rho(0.999)
{
	m_tetMesh = tetMesh;
	m_IsoMaterial = isoMaterial;
	m_vel_old = m_tetMesh->getVelocities();

	n_nodes = m_tetMesh->getNodesNum();
	n_tets = m_tetMesh->getTetsNum();
	m_posk_last = m_posk = m_tetMesh->getNodes();

}


DescentOptimize::~DescentOptimize()
{
}

void DescentOptimize::initialization()
{
	//m_posk = m_pos1 + m_timeStep * m_vel1 + 0.2 * m_timeStep * (m_vel1 - m_vel0);
	m_pos = m_tetMesh->getNodes() + m_h * m_tetMesh->getVelocities();
	m_posk = m_pos + 0.2 * (m_tetMesh->getVelocities()-m_vel_old);
}

Eigen::MatrixXf DescentOptimize::computeGradient()
{
	Eigen::MatrixXf g = m_posk - m_pos;
	for (int i = 0; i < n_nodes; ++i)
	{
		g.col(i) = g.col(i) * m_tetMesh->getMasses()[i] / (m_h*m_h);
	}

	Eigen::MatrixXf f = m_IsoMaterial->computeInnerForceFromPos(m_posk);
	g = g - f;

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
	for (int i = 0; i < n_nodes; ++i)
	{
		E += 0.5 * m_tetMesh->getMasses()[i] * ( p.col(i).squaredNorm() + 9.8 * pos(1,i));
	}
	return E;
}

void DescentOptimize::doDescentOpt()
{
	Eigen::MatrixXf delta_q;
	Eigen::VectorXf H_q;
	float m_h2_inv;
	m_alpha = 0.3;
	m_IsoMaterial->computeElasticEnergyFromPos(m_posk);
	Eigen::VectorXf tmp = m_IsoMaterial->getElasticEnergys();
	initialization();
	delta_q = computeGradient();


	for (int k = 0; k < m_iterations; ++k)
	{
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

		for (int i = 0; i < n_nodes; ++i)
		{
			delta_q(0, i) /= -H_q[3 * i + 0];
			delta_q(1, i) /= -H_q[3 * i + 1];
			delta_q(2, i) /= -H_q[3 * i + 2];
		}

		m_posk_next = m_posk + m_alpha * delta_q;

		if (k % 8 == 0)
		{
			if (k == 0)
			{
				m_energy = m_energy_old = computeTotalEnergy(m_posk_next);
			}
			else
			{
				m_energy = computeTotalEnergy(m_posk_next);
				if (m_energy > m_energy_old)
				{
					m_alpha *= 0.7;
					k -= 8;
					m_iterations -= 8;
					m_posk = m_posk_old;
					continue;
				}
				else
				{
					m_energy_old = m_energy;
					m_posk_old   = m_posk;
				}
			}
		}

		if (k < 10) m_omega = 1.0;
		else if(k == 10) m_omega = 2 / (2 - m_rho*m_rho);
		else m_omega = 4 / (4 - m_rho*m_rho*m_omega);


		m_posk_next = 0.9 * (m_posk_next - m_posk) + m_posk;
		m_posk_next = m_omega * (m_posk_next - m_posk_last) + m_posk_last;

		m_posk_last = m_posk;
		m_posk = m_posk_next;
	}
	m_vel_old = m_tetMesh->getVelocities();
	m_tetMesh->getVelocities() = (m_posk - m_tetMesh->getNodes()) / m_h;
	m_tetMesh->getNodes() = m_posk;
}
