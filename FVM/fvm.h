#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include <qtimer.h>
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

public slots:
	void DoOneStep();
	void DoRun();
	void DoPause();
	void DoLoadConfig();
	void DoStop();

private:
	void initSignalSlotConnections();

	Ui::FVMClass ui;
	std::shared_ptr<TetMesh> m_tetMesh;
	std::shared_ptr<TimeIntegration> m_integrator;

	QTimer m_idleTimer;
};

#endif // FVM_H
