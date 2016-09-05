#pragma once

#include "QGLViewer\qglviewer.h"

class RenderWidget :
	public QGLViewer
{
public:
	RenderWidget(QWidget *parent);
	virtual ~RenderWidget();
protected:
	virtual void draw();
	virtual void init();
};

