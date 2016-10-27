#pragma once

#include "QGLViewer\qglviewer.h"
#include "qopenglfunctions_4_5_core.h"
#include "qopenglcontext.h"
#include <QTime>
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
	void renderText3D(double x, double y, double z, QString text, QPainter *painter);
	void renderText2D(double x, double y, QString text, QPainter *painter);
	void restartTime() { m_time.restart(); }
	void startTime() { m_time.start(); };

protected:
	virtual void draw();
	virtual void init();

	//virtual void paintGL() { update(); }
	virtual void paintEvent(QPaintEvent *event);

private:
	void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);
	GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
		const GLdouble model[16], const GLdouble proj[16],
		const GLint viewport[4],
		GLdouble * winx, GLdouble * winy, GLdouble * winz);

	void drawTestCube();

	QOpenGLFunctions_4_5_Core *render;    // reserve for modern glsl rendering
	std::shared_ptr<TetMesh> gl_tetmesh;

	QTime m_time;
	int m_fps;
};

inline void RenderWidget::transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4])
{
#define M(row,col)  m[col*4+row]
	out[0] =
		M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
	out[1] =
		M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
	out[2] =
		M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
	out[3] =
		M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
#undef M
}

inline GLint RenderWidget::project(GLdouble objx, GLdouble objy, GLdouble objz,
	const GLdouble model[16], const GLdouble proj[16],
	const GLint viewport[4],
	GLdouble * winx, GLdouble * winy, GLdouble * winz)
{
	GLdouble in[4], out[4];

	in[0] = objx;
	in[1] = objy;
	in[2] = objz;
	in[3] = 1.0;
	transformPoint(out, model, in);
	transformPoint(in, proj, out);

	if (in[3] == 0.0)
		return GL_FALSE;

	in[0] /= in[3];
	in[1] /= in[3];
	in[2] /= in[3];

	*winx = viewport[0] + (1 + in[0]) * viewport[2] / 2;
	*winy = viewport[1] + (1 + in[1]) * viewport[3] / 2;

	*winz = (1 + in[2]) / 2;
	return GL_TRUE;
}