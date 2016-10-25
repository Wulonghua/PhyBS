#include "RenderWidget.h"


RenderWidget::RenderWidget(QWidget *parent)
{

}


RenderWidget::~RenderWidget()
{

}

void RenderWidget::init()
{
	restoreStateFromFile();
	this->setSceneRadius(2);

	render = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_4_5_Core>();
	if (!render) {
		std::cerr << "Could not obtain required OpenGL context version";
		exit(1);
	}

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_LIGHTING);
	glDisable(GL_POINT_SMOOTH);


	// Light setup
	//glDisable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	//Light default parameters
	const GLfloat light_ambient[4] = { 0.3, 0.3, 0.3, 1.0 };
	const GLfloat light_diffuse[4] = { 0.3, 0.3, 0.3, 1.0 };
	const GLfloat light_specular[4] = { 0.3, 0.3, 0.3, 1.0 };

	glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.1f);
	glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.3f);
	glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.3f);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
}

void RenderWidget::draw()
{
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glDisable(GL_LIGHT1);
	//drawTestCube();
	gl_tetmesh->drawTetBoundFace();
}

void RenderWidget::postDraw()
{
	renderText2D(10, 30, "FPS: ");
}

void RenderWidget::renderText2D(double x, double y, QString text)
{
	QPainter painter(this);
	painter.setPen(Qt::black);
	painter.setFont(QFont("Helvetica", 20));
	painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	painter.drawText(x, y, text); // z = pointT4.z + distOverOp / 4
	painter.end();
}

void RenderWidget::renderText3D(double x, double y, double z, QString text)
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

	QPainter painter(this);
	painter.setPen(Qt::black);
	painter.setFont(QFont("Helvetica", 10));
	painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
	painter.drawText(textPosX, textPosY, text); // z = pointT4.z + distOverOp / 4
	painter.end();
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