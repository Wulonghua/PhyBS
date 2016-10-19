#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include <qtimer.h>
#include <qlabel.h>
#include <qfiledialog.h>
#include <memory>

#include "ui_fvm.h"
#include "TetMesh.h"
#include "TimeIntegration.h"
#include "IsotropicMaterial.h"
#include "IsotropicNeohookeanMaterial.h"

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
	void DoImportNodes();   // import nodes' positions.
	void DoStop();
	void DoTest();

private:
	void initSignalSlotConnections();

	Ui::FVMClass ui;
	std::shared_ptr<TetMesh> m_tetMesh;
	std::shared_ptr<TimeIntegration> m_integrator;
	std::shared_ptr<IsotropicNeohookeanMaterial> m_IsoMaterial;

	QTimer *m_idleTimer;
	QLabel *m_statusLabel;

	int m_iter;
	int m_frameID;
};

#endif // FVM_H
