#include "fvm.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	//QSurfaceFormat format;
	//format.setRenderableType(QSurfaceFormat::OpenGL);
	//format.setProfile(QSurfaceFormat::CoreProfile);
	//format.setVersion(4, 5);

	FVM w;
	w.show();
	return a.exec();
}
