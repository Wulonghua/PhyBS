#ifndef FVM_H
#define FVM_H

#include <QtWidgets/QMainWindow>
#include "ui_fvm.h"

class FVM : public QMainWindow
{
	Q_OBJECT

public:
	FVM(QWidget *parent = 0);
	~FVM();

private:
	Ui::FVMClass ui;
};

#endif // FVM_H
