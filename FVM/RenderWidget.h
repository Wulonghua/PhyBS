#pragma once

#include "QGLViewer\qglviewer.h"
#include "qopenglfunctions_4_5_core.h"
#include "qopenglcontext.h"
#include "TetMesh.h"

#include <iostream>

class RenderWidget :
	public QGLViewer
{
public:
	RenderWidget(QWidget *parent);
	virtual ~RenderWidget();
protected:
	virtual void draw();
	virtual void init();
private:

	void DrawTestCube();

	QOpenGLFunctions_4_5_Core *render;    // reserve for modern glsl rendering

	std::shared_ptr<TetMesh> gl_tetmesh;
};

