#pragma once

#include <Eigen/Dense>
#include <iostream>

class PosBaseDynamic
{
public:
	PosBaseDynamic(Eigen::MatrixXf& X, int numtets);
	~PosBaseDynamic();

	// [Muller et al. 14] Strain Based Dynamics
	void doStepStrainConstraints(Eigen::MatrixXf &pos,
								 Eigen::MatrixXf &vel,
								 const Eigen::MatrixXf &forces,
								 const Eigen::MatrixXf &invRestMats,
								 const Eigen::MatrixXi &tets,
								 const Eigen::VectorXf &invMass,
								 float t);

private:
	// stretchK and shearK are stiffness coefficients in streatch and shear parts respectively
	void solveStrainConstraints(const Eigen::Vector3f &p0, const float &invMass0, Eigen::Vector3f &d_p0,
								const Eigen::Vector3f &p1, const float &invMass1, Eigen::Vector3f &d_p1,
								const Eigen::Vector3f &p2, const float &invMass2, Eigen::Vector3f &d_p2, 
								const Eigen::Vector3f &p3, const float &invMass3, Eigen::Vector3f &d_p3,
								const Eigen::Matrix3f &invRestMat, const Eigen::Vector3f &stretchK, const Eigen::Vector3f &shearK);

	Eigen::MatrixXf m_pos;	// nodes' positions during iteration
	int n_tets;
	int n_maxIters;
};

