#include "fvm.h"

FVM::FVM(QWidget *parent)
	: QMainWindow(parent), m_iter(0), m_frameID(0), m_numThreads(8)
{
	ui.setupUi(this);
	m_tetMesh = std::make_shared<TetMesh>();
	//m_integrator = std::make_shared<TimeIntegration>(m_tetMesh->getNodesNum(),m_tetMesh->getMasses());
	m_integrator = std::make_shared<TimeIntegration>(m_tetMesh->getNodesNum(),
		m_tetMesh->getMasses(),m_tetMesh->getConstraintIDs(),m_tetMesh->getRestPosition());
	ui.glWidget->setGLTetMesh(m_tetMesh);
	ui.glWidget->startTime();
	m_IsoMaterial = std::make_shared<IsotropicNeohookeanMaterial>(m_tetMesh);


	//int num_nodes, const float *nodes, const float *restPoses, const int *constraintsMask,
	//int num_tets, const int *tets, const float youngs, const float nu, const float density,
	//int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval
	std::vector<int> mask;
	mask.resize(m_tetMesh->getNodesNum());
	std::fill(mask.begin(), mask.end(), 0);
	auto c = m_tetMesh->getConstraintIDs();
	for (int i = 0; i < c.size(); i++)
	{
		mask[c[i]] = 1;
	}
	auto K = m_IsoMaterial->getGlobalK();
	m_cudaInterface = std::make_shared<CUDAInterface>(
	m_tetMesh->getNodesNum(), m_tetMesh->getNodes().data(), m_tetMesh->getRestPosition().data(), mask.data(),
	m_tetMesh->getTetsNum(), m_tetMesh->getTets().data(), m_tetMesh->getE(), m_tetMesh->getNu(),1000,
	K.nonZeros(), K.rows(),K.valuePtr(),K.outerIndexPtr(),K.innerIndexPtr(),m_IsoMaterial->getDiagonalIdx().data(),m_IsoMaterial->getKIDinCSRval().data());

	std::cout << "cuda initialized" << std::endl;

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
	connect(ui.actionImportNodes, SIGNAL(triggered()), this, SLOT(DoImportNodes()));
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

void FVM::DoImportNodes()
{
	QFileDialog *fileDialog = new QFileDialog(this);
	fileDialog->setWindowTitle(QStringLiteral("Import Nodes Positions"));
	fileDialog->setNameFilter(QString("Node(*.node)"));

	if (fileDialog->exec() == QDialog::Accepted)
	{
		QString nodeFile = fileDialog->selectedFiles()[0];
		if (nodeFile.isEmpty())
		{
			std::cout << "Failed to update nodes' positions." << std::endl;
			return;
		}
		m_tetMesh->updateNodesFromFile(nodeFile);
		ui.glWidget->update();
		std::cout << "Nodes' positions updated.";
	}
	else
	{
		std::cout << "Failed to update nodes' positions." << std::endl;
	}

}

void FVM::DoOneStep()
{
	ui.glWidget->restartTime();
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
	Eigen::MatrixXf forces = m_IsoMaterial->computeInnerForcesfromFhats2(m_numThreads);
	//Eigen::MatrixXf K = m_IsoMaterial->computeStiffnessMatrix(0);
	Eigen::SparseMatrix<float, Eigen::RowMajor> sK = m_IsoMaterial->computeGlobalStiffnessMatrix(m_numThreads);
	//std::cout << "before integration: " << std::endl;
	//std::cout << "positions: " << std::endl;
	//std::cout << m_tetMesh->getNodes() << std::endl;

	//std::cout << "forces: " << std::endl;
	//std::cout << m_tetMesh->getForces() << std::endl;

	//std::cout << "velocity: " << std::endl;
	//std::cout << m_tetMesh->getVelocities() << std::endl;

	//std::cout << std::endl;

	//m_integrator->BackEuler(m_tetMesh->getNodes(),
	//	m_tetMesh->getVelocities(),
	//	forces, sK);

	m_integrator->BackEuler(m_tetMesh->getNodes(), m_tetMesh->getRestPosition(),
		m_tetMesh->getVelocities(),
		forces, sK);

	//m_tetMesh->getNodes().col(0) = Eigen::Vector3f(0.1, 0.1, 0.1);
	//m_tetMesh->getVelocities().col(0) = Eigen::Vector3f::Zero();
	
	

	//m_nodes.col(1) = Eigen::Vector3f(-1.0, 1.0, -1.0);
	//m_velocities.col(1) = Eigen::Vector3f::Zero();

	//std::cout << "after integration: " << std::endl;
	//std::cout << "positions: " << std::endl;
	//std::cout << m_tetMesh->getNodes() << std::endl;

	//std::cout << "forces: " << std::endl;
	//std::cout << m_tetMesh->getForces() << std::endl;

	//std::cout << "velocity: " << std::endl;
	//std::cout << m_tetMesh->getVelocities() << std::endl;

	//QString node_file = QStringLiteral("bar_%1").arg(++m_frameID) + QStringLiteral(".node");
	//m_tetMesh->writeNodes(node_file);
	//ui.glWidget->saveSnapshot(snap_file, true);
		//QString snap_file = QStringLiteral("Teran_%1").arg(++m_frameID) + QStringLiteral(".bmp");
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
	//std::cout<< m_IsoMaterial->computeInnerForcesfromFhats2();

	//m_tetMesh->writeMatrix("Us_CPU.csv", m_IsoMaterial->getUs());
	//m_tetMesh->writeMatrix("Vs_CPU.csv", m_IsoMaterial->getVs());
	//m_tetMesh->writeMatrix("Fhats.csv", m_IsoMaterial->getFhats());
	//Eigen::SparseMatrix<float, Eigen::RowMajor> sK = m_IsoMaterial->computeGlobalStiffnessMatrix();
	//Eigen::MatrixXf K(sK);
	//m_tetMesh->writeMatrix("K.csv", K);
	//auto m = m_IsoMaterial->computeStiffnessMatrix(0);
	//std::cout << m;
	//Eigen::MatrixXf forces = m_IsoMaterial->computeInnerForcesfromFhats();
	//Eigen::MatrixXf K = m_IsoMaterial->computeStiffnessMatrix(0);
	//std::cout << "K: " << std::endl;
	//std::cout << K << std::endl;

	/****************** test for cpu parallel elapse time *************************
	int elapse;

	ui.glWidget->restartTime();
	Eigen::MatrixXf forces = m_IsoMaterial->computeInnerForcesfromFhats(m_numThreads);

	elapse = ui.glWidget->restartTime();
	std::cout << "time to compute force: " << elapse << std::endl;


	Eigen::SparseMatrix<float, Eigen::RowMajor> sK = m_IsoMaterial->computeGlobalStiffnessMatrix(m_numThreads);

	elapse = ui.glWidget->restartTime();
	std::cout << "time to compute K: " << elapse << std::endl;

	m_integrator->BackEuler(m_tetMesh->getNodes(), m_tetMesh->getRestPosition(),
		m_tetMesh->getVelocities(),
		forces, sK);
	elapse = ui.glWidget->restartTime();
	std::cout << "time to solve the system: " << elapse << std::endl;
	ui.glWidget->update();
	**************************************************************************************/
	//int num_nodes, const float *nodes, const float *restPoses, const int *constraintsMask,
	//int num_tets, const int *tets, const float youngs, const float nu, const float density,
	//int csr_nnz, int csr_m, float *csr_val, int *csr_row, int *csr_col, int *csr_diagonalIdx, int *csr_kIDinCSRval
	//std::vector<int> mask;
	//mask.resize(m_tetMesh->getNodesNum());
	//std::fill(mask.begin(), mask.end(), 0);
	//auto c = m_tetMesh->getConstraintIDs();
	//for (int i = 0; i < c.size(); i++)
	//{
	//	mask[c[i]] = 1;
	//}
	//auto K = m_IsoMaterial->getGlobalK();
	//m_cudaInterface = std::make_shared<CUDAInterface>(
	//m_tetMesh->getNodesNum(), m_tetMesh->getNodes().data(), m_tetMesh->getRestPosition().data(), mask.data(),
	//m_tetMesh->getTetsNum(), m_tetMesh->getTets().data(), m_tetMesh->getE(), m_tetMesh->getNu(),1000,
	//K.nonZeros(), K.rows(),K.valuePtr(),K.outerIndexPtr(),K.innerIndexPtr(),m_IsoMaterial->getDiagonalIdx().data(),m_IsoMaterial->getKIDinCSRval().data());

	//std::cout << "cuda initialized" << std::endl;
	m_cudaInterface->computeInnerforces();

	m_cudaInterface->computeGlobalStiffnessMatrix();


}
