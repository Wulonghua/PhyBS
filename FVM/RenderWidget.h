#pragma once

#include "QGLViewer\qglviewer.h"
#include "qopenglfunctions_4_5_core.h"
#include "qopenglcontext.h"
#include "qopenglbuffer.h"
#include "qopenglshaderprogram.h"
#include "qdebug.h"
#include "TetMesh.h"

#include <iostream>
#include <memory>

class RenderWidget :
	public QGLViewer
{
public:
	RenderWidget(QWidget *parent);
	virtual ~RenderWidget();

	void setGLTetMesh(std::shared_ptr<TetMesh> tet) { gl_tetmesh = tet; }

protected:
	virtual void draw();
	virtual void init();

private:

	void drawTestCube();
	void initTestShaders();
	void drawTestShaders();
	void printVersionInformation();


	QOpenGLFunctions_4_5_Core *render;    // reserve for modern glsl rendering
	std::shared_ptr<TetMesh> gl_tetmesh;

	QOpenGLShaderProgram m_tetShaderProgram;
	QOpenGLBuffer m_NodePosBuf;
	QOpenGLBuffer m_NodeIdxBuf;
	QOpenGLBuffer m_NodeColorBuf;

};

