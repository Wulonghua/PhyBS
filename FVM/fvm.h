#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include "ui_fvm.h"
#include "TetMesh.h"

#include <memory>

class FVM : public QMainWindow
{
	Q_OBJECT

public:
	FVM(QWidget *parent = 0);
	~FVM();

private:
	Ui::FVMClass ui;
	std::shared_ptr<TetMesh> m_tetMesh;


};

#endif // FVM_H
