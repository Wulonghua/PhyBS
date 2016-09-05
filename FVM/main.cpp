#include "fvm.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	FVM w;
	w.show();
	return a.exec();
}
