#include "fvm.h"

FVM::FVM(QWidget *parent)
	: QMainWindow(parent)
{
	m_tetMesh = std::make_shared<TetMesh>();
	ui.setupUi(this);
	ui.glWidget->setGLTetMesh(m_tetMesh);
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
}

void FVM::DoLoadConfig()
{

}

void FVM::DoOneStep()
{

	ui.glWidget->update();
}

void FVM::DoRun()
{

}

void FVM::DoPause()
{

}

void FVM::DoStop()
{

}
