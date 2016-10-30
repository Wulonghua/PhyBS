#include "RenderWidget.h"
#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

RenderWidget::RenderWidget(QWidget *parent) : m_fps(0), m_elapses(0), m_iter(0), m_iterMax(10), m_picked(false)
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
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	gl_tetmesh->drawTetBoundFace();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
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
		qglviewer::Vec orig;
		QPoint p;
		Eigen::Vector3d origin, direction, endp;

		//std::cout << "width: " << this->width() << std::endl;
		//std::cout << "height: " << this->height() << std::endl;
		//std::cout << "x: " << e->localPos().x() << std::endl;
		//std::cout << "y: " << e->localPos().y() << std::endl;
	
		double matModelView[16], matProjection[16];
		int viewport[4];

		// get matrix and viewport:
		glGetDoublev(GL_MODELVIEW_MATRIX, matModelView);
		glGetDoublev(GL_PROJECTION_MATRIX, matProjection);
		glGetIntegerv(GL_VIEWPORT, viewport);

		// window pos of mouse, Y is inverted on Windows
		double winX = e->pos().x();
		double winY = viewport[3] - e->pos().y();

		// get point on the 'far' plane (third param is set to 1.0)
		gluUnProject(winX, winY, 1.0, matModelView, matProjection,
			viewport, &endp[0], &endp[1], &endp[2]);

		orig = camera()->position();
		
		origin << orig.x, orig.y, orig.z;
		direction = (endp - origin).normalized();

		std::cout << "picked: "<<gl_tetmesh->pickFacebyRay(origin, direction)<<std::endl;

		
		std::cout << "orig: " << orig.x << ","<< orig.y<< "," << orig.z << std::endl;
		std::cout << "dir: " << direction[0] << "," << direction[1] << "," << direction[2] << std::endl;
		

	}
	else
	{
		QGLViewer::mousePressEvent(e);
	}


}

void RenderWidget::renderText2D(double x, double y, QString text, QPainter *painter)
{
	painter->setPen(Qt::darkCyan);
	painter->setFont(QFont("Arial", 20));
	painter->setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	painter->drawText(x, y, text); // z = pointT4.z + distOverOp / 4
}

void RenderWidget::renderText3D(double x, double y, double z, QString text, QPainter *painter)
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