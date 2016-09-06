#include "fvm.h"

FVM::FVM(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	m_tetMesh = std::make_shared<TetMesh>();
}

FVM::~FVM()
{

}
