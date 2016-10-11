#include "fvm.h"

FVM::FVM(QWidget *parent)
	: QMainWindow(parent), m_iter(0), m_frameID(0)
{
	ui.setupUi(this);
	m_tetMesh = std::make_shared<TetMesh>();
	m_integrator = std::make_shared<TimeIntegration>(m_tetMesh->getNodesNum(),m_tetMesh->getMasses());
	ui.glWidget->setGLTetMesh(m_tetMesh);
	m_IsoMaterial = std::make_shared<IsotropicNeohookeanMaterial>(m_tetMesh);

	m_statusLabel = new QLabel(this);
	m_statusLabel->setText(QStringLiteral(" #Nodes: %1	#Tets: %2	Time step: %3s")
					.arg(m_tetMesh->getNodesNum())
					.arg(m_tetMesh->getTetsNum())
					.arg(m_integrator->getTimeStep()));
	ui.statusBar->addWidget(m_statusLabel);

	m_idleTimer = new QTimer(this);
	initSignalSlotConnections();
}

FVM::~FVM()
{

}

void FVM::initSignalSlotConnections()
{
	connect(ui.actionLoadConfig, SIGNAL(triggered()), this, SLOT(DoLoadConfig()));
	connect(ui.actionStep, SIGNAL(triggered()), this, SLOT(DoOneStep()));
	connect(ui.actionRun, SIGNAL(triggered()), this, SLOT(DoRun()));
	connect(ui.actionPause, SIGNAL(triggered()), this, SLOT(DoPause()));
	connect(ui.actionStop, SIGNAL(triggered()), this, SLOT(DoStop()));
	connect(ui.actionTest, SIGNAL(triggered()), this, SLOT(DoTest()));
	connect(m_idleTimer,SIGNAL(timeout()),this, SLOT(DoOneStep()));
}

void FVM::DoLoadConfig()
{

}

void FVM::DoOneStep()
{
	/************************linear model using explicit time integration*****************************/
	//if (m_iter++ < 10000)
	//{
	//	m_tetMesh->computeForces();
	//	m_integrator->simuExplicit(m_tetMesh->getNodes(),
	//		m_tetMesh->getVelocities(),
	//		m_tetMesh->getForces(),
	//		m_tetMesh->getMasses());
	//	m_tetMesh->updateNodesVelocities(m_integrator->getPositions(), m_integrator->getVelocities());
	//	ui.glWidget->update();
	//}
	//else
	//{
	//	m_iter = 0;
	//	QString snap_file = QStringLiteral("snapshot_%1").arg(++m_frameID) + QStringLiteral(".png");
	//	ui.glWidget->saveSnapshot(snap_file, true);
	//}
	/**********************************************************************************************/

	/******************************Back Euler integration**********************************/
	Eigen::MatrixXd forces = m_IsoMaterial->computeInnerForcesfromFhats();
	//Eigen::MatrixXd K = m_IsoMaterial->computeStiffnessMatrix(0);
	Eigen::SparseMatrix<double> sK = m_IsoMaterial->computeGlobalStiffnessMatrix();
	//std::cout << "before integration: " << std::endl;
	//std::cout << "positions: " << std::endl;
	//std::cout << m_tetMesh->getNodes() << std::endl;

	//std::cout << "forces: " << std::endl;
	//std::cout << m_tetMesh->getForces() << std::endl;

	//std::cout << "velocity: " << std::endl;
	//std::cout << m_tetMesh->getVelocities() << std::endl;

	//std::cout << std::endl;

	m_integrator->BackEuler(m_tetMesh->getNodes(),
		m_tetMesh->getVelocities(),
		forces, sK);

	//m_tetMesh->getNodes().col(0) = Eigen::Vector3d(0.1, 0.1, 0.1);
	//m_tetMesh->getVelocities().col(0) = Eigen::Vector3d::Zero();
	
	

	//m_nodes.col(1) = Eigen::Vector3d(-1.0, 1.0, -1.0);
	//m_velocities.col(1) = Eigen::Vector3d::Zero();

	//std::cout << "after integration: " << std::endl;
	//std::cout << "positions: " << std::endl;
	//std::cout << m_tetMesh->getNodes() << std::endl;

	//std::cout << "forces: " << std::endl;
	//std::cout << m_tetMesh->getForces() << std::endl;

	//std::cout << "velocity: " << std::endl;
	//std::cout << m_tetMesh->getVelocities() << std::endl;

	//std::cout << std::endl;
	//QString snap_file = QStringLiteral("torus_%1").arg(++m_frameID) + QStringLiteral(".bmp");
	//ui.glWidget->saveSnapshot(snap_file, true);
	/***********************************************************************************************/
	ui.glWidget->update();

}

void FVM::DoRun()
{
	m_idleTimer->start(0);

}

void FVM::DoPause()
{
	m_idleTimer->stop();
}

void FVM::DoStop()
{
	
}

void FVM::DoTest()
{
	//m_IsoMaterial->computeInnerForcesfromFhats();
	//auto m = m_IsoMaterial->computeStiffnessMatrix(0);
	//std::cout << m;
	Eigen::MatrixXd forces = m_IsoMaterial->computeInnerForcesfromFhats();
	Eigen::MatrixXd K = m_IsoMaterial->computeStiffnessMatrix(0);
	std::cout << "K: " << std::endl;
	std::cout << K << std::endl;

	Eigen::SparseMatrix<double> gK = m_IsoMaterial->computeGlobalStiffnessMatrix();
	Eigen::MatrixXd Kt;
	Kt = Eigen::MatrixXd(gK);
	std::cout << "Kt: " << std::endl;
	std::cout << Kt << std::endl;
}
