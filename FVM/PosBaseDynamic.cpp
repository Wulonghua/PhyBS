#include "PosBaseDynamic.h"

PosBaseDynamic::PosBaseDynamic(Eigen::MatrixXf& X, int numtets):
m_pos(X), n_tets(numtets), n_maxIters(5), m_eps(1e-6)
{
}

PosBaseDynamic::~PosBaseDynamic()
{
}

void PosBaseDynamic::doStepStrainConstraints(Eigen::MatrixXf &pos,
										     Eigen::MatrixXf &vel,
										     const Eigen::MatrixXf &forces,
										     const Eigen::MatrixXf &invRestMat,
										     const Eigen::MatrixXi &tets,
										     const Eigen::VectorXf &invMass,
										     float t)
{
	Eigen::Vector3f stretchStiff(0.8,0.8,0.8);
	Eigen::Vector3f shearStiff(0.8,0.8,0.8);
	Eigen::Vector3f d_p0, d_p1, d_p2, d_p3;
	int idx[4]; // node's index
	int n_nodes = invMass.size();
	Eigen::MatrixXf d_v = Eigen::MatrixXf::Zero(3, n_nodes);
	for (int i = 0; i < n_nodes; ++i)
	{
		d_v.col(i) = t * invMass(i) * forces.col(i);
	}
	m_pos = pos + t * (vel + d_v);
	for (int i = 0; i < n_maxIters; ++i)
	{
		for (int j = 0; j < n_tets; ++j)
		{
			idx[0] = tets(0, j);
			idx[1] = tets(1, j);
			idx[2] = tets(2, j);
			idx[3] = tets(3, j);
			solveStrainConstraints(
				m_pos.col(idx[0]), invMass(idx[0]), d_p0,
				m_pos.col(idx[1]), invMass(idx[1]), d_p1,
				m_pos.col(idx[2]), invMass(idx[2]), d_p2,
				m_pos.col(idx[3]), invMass(idx[3]), d_p3,
				invRestMat.block<3, 3>(0, 3 * j), stretchStiff, shearStiff);
			m_pos.col(idx[0]) += d_p0;
			m_pos.col(idx[1]) += d_p1;
			m_pos.col(idx[2]) += d_p2;
			m_pos.col(idx[3]) += d_p3;
		}
	}
	vel = (m_pos - pos) / t;
	pos = m_pos;
}

void PosBaseDynamic::solveStrainConstraints(const Eigen::Vector3f &p0, const float &invMass0, Eigen::Vector3f &d_p0,
											const Eigen::Vector3f &p1, const float &invMass1, Eigen::Vector3f &d_p1,
											const Eigen::Vector3f &p2, const float &invMass2, Eigen::Vector3f &d_p2,
											const Eigen::Vector3f &p3, const float &invMass3, Eigen::Vector3f &d_p3,
											const Eigen::Matrix3f &invRestMat, const Eigen::Vector3f &stretchK, const Eigen::Vector3f &shearK)
{
	d_p0.setZero();
	d_p1.setZero();
	d_p2.setZero();
	d_p3.setZero();

	Eigen::Matrix3f P;
	Eigen::Vector3f fi, fj;
	Eigen::Vector3f d[4];
	float Sij, lambda;


	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			P.col(0) = (p1 + d_p1) - (p0 + d_p0);
			P.col(1) = (p2 + d_p2) - (p0 + d_p0);
			P.col(2) = (p3 + d_p3) - (p0 + d_p0);

			fi = P * invRestMat.col(i);
			fj = P * invRestMat.col(j);
			Sij = fi.dot(fj);

			d[0] = Eigen::Vector3f::Zero();

			for (int k = 0; k < 3; ++k)
			{
				d[k + 1] = fj * invRestMat(k, i) + fi * invRestMat(k, j);
				d[0] -= d[k + 1];
			}

			lambda = invMass0 * d[0].squaredNorm() +
					 invMass1 * d[1].squaredNorm() +
					 invMass2 * d[2].squaredNorm() +
					 invMass3 * d[3].squaredNorm();

			if (lambda < m_eps)
				continue;

			if (i == j)
			{
				lambda = (Sij - 1.0) / lambda * stretchK[i];
			}
			else
			{
				lambda = Sij / lambda * shearK[i + j - 1];
			}

			d_p0 -= lambda * invMass0 * d[0];
			d_p1 -= lambda * invMass1 * d[1];
			d_p2 -= lambda * invMass2 * d[2];
			d_p3 -= lambda * invMass3 * d[3];

		}
	}

}
