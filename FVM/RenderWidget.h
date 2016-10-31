#pragma once

#include "QGLViewer\qglviewer.h"
#include "qopenglfunctions_4_5_core.h"
#include "qopenglcontext.h"
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QKeyEvent>
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
	int restartTime() { return m_time.restart(); }
	void startTime() { m_time.start(); };

protected:

	virtual void draw();
	virtual void init();

	//virtual void paintGL() { update(); }
	virtual void paintEvent(QPaintEvent *e);
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);

private:
	void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);
	GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
		const GLdouble model[16], const GLdouble proj[16],
		const GLint viewport[4],
		GLdouble * winx, GLdouble * winy, GLdouble * winz);

	void drawTestCube();

	QOpenGLFunctions_4_5_Core *render;    // reserve for modern glsl rendering
	std::shared_ptr<TetMesh> gl_tetmesh;

	QElapsedTimer m_time;
	int m_fps;
	int m_elapses;
	int m_iter;
	const int m_iterMax;

	double m_ModelView[16];
	double m_Projection[16];
	int m_viewport[4];
	Eigen::Map<Eigen::Matrix<double,4,4>> m_matModelView;
	Eigen::Map<Eigen::Matrix<double, 4, 4>> m_matProjection;
	//Eigen::Matrix4d m_matModelViewInverse;
	//Eigen::Matrix4d m_matProjectionInverse;

	bool m_picked;
	int m_picki;
	Eigen::Vector3d m_line_end1;
	Eigen::Vector3d m_line_end2;
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