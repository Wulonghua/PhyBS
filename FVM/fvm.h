#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include <memory>

#include "ui_fvm.h"
#include "TetMesh.h"
#include "TimeIntegration.h"

class FVM : public QMainWindow
{
	Q_OBJECT

public:
	FVM(QWidget *parent = 0);
	~FVM();

private:
	Ui::FVMClass ui;
	std::shared_ptr<TetMesh> m_tetMesh;
	std::shared_ptr<TimeIntegration> m_integrator;
};

#endif // FVM_H
