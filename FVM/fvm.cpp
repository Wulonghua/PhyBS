#include "fvm.h"

FVM::FVM(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	m_tetMesh = std::make_shared<TetMesh>();
	m_integrator = std::make_shared<TimeIntegration>(m_tetMesh->getNodesNum());
	ui.glWidget->setGLTetMesh(m_tetMesh);
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
	connect(m_idleTimer,SIGNAL(timeout()),this, SLOT(DoOneStep()));
}

void FVM::DoLoadConfig()
{

}

void FVM::DoOneStep()
{
	m_tetMesh->computeForces();
	m_integrator->simuExplicit(	m_tetMesh->getNodes(),
								m_tetMesh->getVelocities(),
								m_tetMesh->getForces(),
								m_tetMesh->getMasses());
	m_tetMesh->updateNodesVelocities(m_integrator->getPositions(), m_integrator->getVelocities());
	ui.glWidget->update();
}

void FVM::DoRun()
{
	m_idleTimer->start(0);

}

void FVM::DoPause()
{
	
}

void FVM::DoStop()
{
	m_idleTimer->stop();
}
