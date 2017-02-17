#include "DescentOptimize.h"


DescentOptimize::DescentOptimize(std::shared_ptr<TetMesh> tetMesh, std::shared_ptr<IsotropicNeohookeanMaterial> isoMaterial)
{
	m_tetMesh = tetMesh;
	m_IsoMaterial = isoMaterial;
	m_vel_old = m_tetMesh->getVelocities();
}


DescentOptimize::~DescentOptimize()
{
}

void DescentOptimize::initialization()
{
	//m_posk = m_pos1 + m_timeStep * m_vel1 + 0.2 * m_timeStep * (m_vel1 - m_vel0);
	m_posk = m_tetMesh->getNodes() + m_timeStep * m_tetMesh->getVelocities() + 0.2 * (m_tetMesh->getVelocities()-m_vel_old);
}

Eigen::MatrixXf DescentOptimize::computeGradient()
{
	Eigen::MatrixXf g = m_posk - m_tetMesh->getNodes() - m_timeStep*m_tetMesh->getVelocities();

}
