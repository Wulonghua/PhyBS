#include "RenderWidget.h"
#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

RenderWidget::RenderWidget(QWidget *parent) :
m_fps(0), m_elapses(0), m_iter(0), m_iterMax(10), m_picked(false),
m_matModelView(m_ModelView), m_matProjection(m_Projection)
{

}


RenderWidget::~RenderWidget()
{

}

void RenderWidget::init()
{
	restoreStateFromFile();
	this->setSceneRadius(5);

	render = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_4_5_Core>();
	if (!render)
	{
		std::cerr << "Could not obtain required OpenGL context version";
		exit(1);
	}
}

void RenderWidget::draw()
{
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	gl_tetmesh->drawTetBoundFace();
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if (m_picked)
	{
		m_line_end1 = gl_tetmesh->getFaceCenter(m_picki);
		gl_tetmesh->drawDraggedNodes(m_picki);
		glColor3d(1, 0, 0);
		glBegin(GL_LINES);
		glVertex3fv(m_line_end1.data());
		glVertex3fv(m_line_end2.data());
		glEnd();
	}
		
}


void RenderWidget::paintEvent(QPaintEvent *e)
{
	Q_UNUSED(e)
	QPainter painter(this);
	painter.beginNativePainting();
	// Save current OpenGL state
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	// Reset OpenGL parameters
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_MULTISAMPLE);
	static GLfloat lightPosition[4] = { 3.0, 3.0, 3.0, 1.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	glClearColor(1.0, 1.0, 1.0, 1.0);
	//glDisable(GL_POINT_SMOOTH);

	// Light setup
	glEnable(GL_LIGHT1);

	//Light default parameters
	const GLfloat light_ambient[4] = { 0.1, 0.1, 0.1, 1.0 };
	const GLfloat light_diffuse[4] = { 0.1, 0.1, 0.1, 1.0 };
	const GLfloat light_specular[4] = { 0.1, 0.1, 0.1, 1.0 };

	glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.1f);
	glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.3f);
	glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.3f);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);

	// Classical 3D drawing, usually performed by paintGL().
	preDraw();
	draw();
	postDraw();
	// Restore OpenGL state
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();

	painter.endNativePainting();

	int elapse = m_time.restart();

	m_elapses += elapse;
	if (++m_iter >= m_iterMax)
	{
		m_fps = m_elapses > 0 ? 1000 * m_iterMax / m_elapses : 60;
		m_iter = 0;
		m_elapses = 0;
	}
	
	renderText2D(10,30,QStringLiteral("FPS: %1").arg(m_fps),&painter);

	painter.end();
}

void RenderWidget::mousePressEvent(QMouseEvent *e)
{
	
	if ((e->button()==Qt::LeftButton) && (e->modifiers()==Qt::ControlModifier))
	{
		qglviewer::Vec orig,dir;
		Eigen::Vector3d origin, direction, endp1,endp2;

		camera()->getModelViewMatrix(m_ModelView);
		camera()->getProjectionMatrix(m_Projection);
		camera()->getViewport(m_viewport);

		//std::cout << "model view: " << m_matModelView << std::endl;
		//std::cout << "project : " << m_matProjection << std::endl;

		//std::cout << "width and height: "<< width() << " " << height()<<std::endl;
		//std::cout << "point: " << e->pos().x() << " " << e->pos().y() << std::endl;

		// window pos of mouse, Y is inverted on Windows
		float winX = e->pos().x();
		float winY = e->pos().y();

		// get point on the 'near' plane 
		gluUnProject(winX, winY, 0.0, m_ModelView, m_Projection,
		             m_viewport, &endp1[0], &endp1[1], &endp1[2]);

		gluUnProject(winX, winY, 1.0, m_ModelView, m_Projection, m_viewport, &endp2[0], &endp2[1], &endp2[2]);
	
		orig = camera()->position();
		//camera()->convertClickToLine(e->pos(), orig, dir);
		origin << orig.x, orig.y, orig.z;
		//direction << dir.x, dir.y, dir.z;
		direction = (endp2 - endp1).normalized();
		m_picki = gl_tetmesh->pickFacebyRay(origin.cast<float>(), direction.cast<float>());
		std::cout << "picked: "<< m_picki <<std::endl;

		if (m_picki > -1)
		{
			//m_line_end1 = origin;
			//m_line_end2 = origin + 20 * direction;
			m_line_end1 = gl_tetmesh->getFaceCenter(m_picki);
			m_line_end2 = m_line_end1;
			m_picked = true;
		}

		//std::cout << "orig: " << orig.x << ","<< orig.y<< "," << orig.z << std::endl;
		//std::cout << "dir: " << direction[0] << "," << direction[1] << "," << direction[2] << std::endl;
	}
	else
	{
		QGLViewer::mousePressEvent(e);
	}
}

void RenderWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (m_picked)
	{
		Eigen::Vector3f A2= gl_tetmesh->getFaceCenter(m_picki);
		m_line_end1 = A2;
		qglviewer::Vec qA2(A2[0], A2[1], A2[2]);
		//qglviewer::Vec qA1 = camera()->cameraCoordinatesOf(qA2);
		qglviewer::Vec qa = camera()->projectedCoordinatesOf(qA2);
		qglviewer::Vec qb;
		qb.x = e->pos().x();
		qb.y = e->pos().y();
		qb.z = qa.z;
		qglviewer::Vec qB2 = camera()->unprojectedCoordinatesOf(qb);
		m_line_end2[0] = qB2.x;
		m_line_end2[1] = qB2.y;
		m_line_end2[2] = qB2.z;

		//gl_tetmesh->dragFace(m_picki, m_line_end2 - m_line_end1);
		gl_tetmesh->dragFaceRing(m_picki, m_line_end2 - m_line_end1);
	}
	else
	{
		QGLViewer::mouseMoveEvent(e);
	}
	update();
}

void RenderWidget::mouseReleaseEvent(QMouseEvent *e)
{
	if (m_picked)
	{
		m_picked = false;
		gl_tetmesh->resetExternalForce();
	}
	QGLViewer::mouseReleaseEvent(e);
}

void RenderWidget::renderText2D(float x, float y, QString text, QPainter *painter)
{
	painter->setPen(Qt::darkCyan);
	painter->setFont(QFont("Arial", 20));
	painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	painter->drawText(x, y, text); // z = pointT4.z + distOverOp / 4
}

void RenderWidget::renderText3D(float x, float y, float z, QString text, QPainter *painter)
{
	int width = this->width();
	int height = this->height();

	GLdouble model[4][4], proj[4][4];
	GLint view[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, &model[0][0]);
	glGetDoublev(GL_PROJECTION_MATRIX, &proj[0][0]);
	glGetIntegerv(GL_VIEWPORT, &view[0]);
	GLdouble textPosX = 0, textPosY = 0, textPosZ = 0;

	project(x, y, z, &model[0][0], &proj[0][0], &view[0], &textPosX, &textPosY, &textPosZ);

	textPosY = height - textPosY; // y is inverted

	painter->setPen(Qt::black);
	painter->setFont(QFont("Helvetica", 10));
	painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	painter->drawText(textPosX, textPosY, text); // z = pointT4.z + distOverOp / 4
}

void RenderWidget::drawTestCube()
{
	// Place light at camera position
	const qglviewer::Vec cameraPos = camera()->position();
	//std::cout<<cameraPos[0]<<" "<<cameraPos[1]<<" "<<cameraPos[2]<<std::endl;
	const GLfloat pos[4] = { cameraPos[0], cameraPos[1], cameraPos[2], 1.0 };
	glLightfv(GL_LIGHT1, GL_POSITION, pos);

	float vertices[8][3] = { { -0.500000, -0.500000, -0.500000 },
	{ 0.500000, -0.500000, -0.500000 },
	{ -0.500000, 0.500000, -0.500000 },
	{ 0.500000, 0.500000, -0.500000 },
	{ -0.500000, -0.500000, 0.500000 },
	{ 0.500000, -0.500000, 0.500000 },
	{ -0.500000, 0.500000, 0.500000 },
	{ 0.500000, 0.500000, 0.500000 } };

	int indices[12][3] = { { 2, 1, 0 }, { 1, 2, 3 },
	{ 4, 2, 0 }, { 2, 4, 6 },
	{ 1, 4, 0 }, { 4, 1, 5 },
	{ 6, 5, 7 }, { 5, 6, 4 },
	{ 3, 6, 7 }, { 6, 3, 2 },
	{ 5, 3, 7 }, { 3, 5, 1 } };

	float normals[12][3] = { { 0, 0, -1 }, {0,0,-1},
							 { -1, 0, 0 }, { -1, 0, 0 },
							 { 0, -1, 0 }, { 0, -1, 0 },
							 { 0, 0, 1 }, { 0, 0, 1 },
							 { 0, 1, 0 }, { 0, 1, 0 },
							 { 1, 0, 0 }, {1,0,0} };

	glColor3f(0.251, 0.424, 0.7);
	for (int i = 0; i < 12; ++i)
	{
		glBegin(GL_TRIANGLES);
		glNormal3fv(&normals[i][0]);
		glVertex3fv(&vertices[indices[i][0]][0]);
		glVertex3fv(&vertices[indices[i][1]][0]);
		glVertex3fv(&vertices[indices[i][2]][0]);
		glEnd();
	}
	
}