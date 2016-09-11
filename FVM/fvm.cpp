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
