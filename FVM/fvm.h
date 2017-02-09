#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include <qtimer.h>
#include <qlabel.h>
#include <qfiledialog.h>
#include <qcombobox.h>
#include <memory>

#include "ui_fvm.h"
#include "TetMesh.h"
#include "TimeIntegration.h"
#include "IsotropicMaterial.h"
#include "IsotropicNeohookeanMaterial.h"
#include "PosBaseDynamic.h"
#include "ProjDynamic.h"
#include "CUDAInterface.h"

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
	void DoExportNodes();
	void DoStop();
	void DoTest();
	void SetTypeComputing(int type);

private:
	void initSignalSlotConnections();

	Ui::FVMClass ui;
	std::shared_ptr<TetMesh> m_tetMesh;
	std::shared_ptr<TimeIntegration> m_integrator;
	std::shared_ptr<IsotropicNeohookeanMaterial> m_IsoMaterial;
	std::shared_ptr<PosBaseDynamic> m_pbd;
	std::shared_ptr<CUDAInterface> m_cudaInterface;

	QTimer *m_idleTimer;
	QLabel *m_statusLabel;

	int m_iter;
	int m_frameID;

	int m_numThreads;

	int m_typeComputing;  // 0: using CPU_force; 1: using GPU_force;

	// UI
	QComboBox *m_comboboxType;

};

#endif // FVM_H
