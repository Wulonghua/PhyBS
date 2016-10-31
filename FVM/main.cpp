#include "fvm.h"
#include <QtWidgets/QApplication>
#include <Eigen/Core>

int main(int argc, char *argv[])
{
	//Eigen::initParallel();
	QApplication a(argc, argv);
	FVM w;
	w.show();
	return a.exec();
}
