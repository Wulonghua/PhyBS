#include "RenderWidget.h"


RenderWidget::RenderWidget(QWidget *parent)
{

}


RenderWidget::~RenderWidget()
{
	m_NodeIdxBuf.destroy();
	m_NodePosBuf.destroy();
}

void RenderWidget::init()
{
	restoreStateFromFile();
	this->setSceneRadius(2);

	//render = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_4_5_Core>();
	render = this->context()->versionFunctions<QOpenGLFunctions_4_5_Core>();
	if (!render) {
		std::cerr << "Could not obtain required OpenGL context version";
		exit(1);
	}

	render->initializeOpenGLFunctions();
	initTestShaders();
	//printVersionInformation();

	render->glClearColor(1.0, 1.0, 1.0, 1.0);

	//glClearColor(1.0, 1.0, 1.0, 1.0);
	//glEnable(GL_LIGHTING);
	//glDisable(GL_POINT_SMOOTH);


	//// Light setup
	////glDisable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);

	////Light default parameters
	//const GLfloat light_ambient[4] = { 1.0, 1.0, 1.0, 1.0 };
	//const GLfloat light_diffuse[4] = { 1.0, 1.0, 1.0, 1.0 };
	//const GLfloat light_specular[4] = { 0.3, 0.3, 0.3, 1.0 };

	//glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION, 0.1f);
	//glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.3f);
	//glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.3f);
	//glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
}

void RenderWidget::draw()
{
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glDisable(GL_LIGHT1);
	//drawTestCube();
	//gl_tetmesh->drawTetBoundFace();

	drawTestShaders();
}

void RenderWidget::initTestShaders()
{	
	m_NodeIdxBuf.create();
	m_NodePosBuf.create();
	m_NodeColorBuf.create();

	GLdouble vertices[] = { -1.0f, -1.0f, 0.0f,
							1.0f, -1.0f, 0.0f,
							0.0f, 1.0f, 0.0f };
	GLdouble colors[] = { 1.0f, 0.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f, 1.0f};

	GLuint indices[] = { 0, 1, 2 };

	m_tetShaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, ".//shaders//test.vert");
	m_tetShaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, ".//shaders//test.frag");
	m_tetShaderProgram.link();
	m_tetShaderProgram.bind();

	m_NodePosBuf.bind();
	m_NodePosBuf.allocate(vertices, 9 * sizeof(GLdouble));

	m_NodeIdxBuf.bind();
	m_NodeIdxBuf.allocate(indices, 3 * sizeof(GLuint));

	m_NodeColorBuf.bind();
	m_NodeColorBuf.allocate(colors, 12 * sizeof(GLdouble));
}

void RenderWidget::drawTestShaders()
{
	m_tetShaderProgram.bind();
	float mvp[16];
	camera()->getModelViewProjectionMatrix(mvp);
	QMatrix4x4 mvp_mat(mvp);
	m_tetShaderProgram.setUniformValue("mvpMatrix", mvp_mat.transposed());

	m_NodeColorBuf.bind();
	int nodeColor = m_tetShaderProgram.attributeLocation("color");
	m_tetShaderProgram.enableAttributeArray(nodeColor);
	m_tetShaderProgram.setAttributeBuffer(nodeColor, GL_DOUBLE, 0, 4, 0);

	m_NodeIdxBuf.bind();
	m_NodePosBuf.bind();
	int nodeLocation = m_tetShaderProgram.attributeLocation("position");
	m_tetShaderProgram.enableAttributeArray(nodeLocation);
	m_tetShaderProgram.setAttributeBuffer(nodeLocation, GL_DOUBLE, 0, 3, 0);


	//GLuint indices[] = { 0, 1, 2 };
	m_NodeIdxBuf.bind();
	render->glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (void*)0);
	//render->glDrawArrays(GL_TRIANGLES, 0, 3);
	m_tetShaderProgram.release();
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

void RenderWidget::printVersionInformation()
{
	QString glType;
	QString glVersion;
	QString glProfile;

	// Get Version Information
	glType = (context()->isOpenGLES()) ? "OpenGL ES" : "OpenGL";
	glVersion = reinterpret_cast<const char*>(glGetString(GL_VERSION));

	// Get Profile Information
#define CASE(c) case QSurfaceFormat::c: glProfile = #c; break
	switch (format().profile())
	{
		CASE(NoProfile);
		CASE(CoreProfile);
		CASE(CompatibilityProfile);
	}
#undef CASE

	// qPrintable() will print our QString w/o quotes around it.
	qDebug() << qPrintable(glType) << qPrintable(glVersion) << "(" << qPrintable(glProfile) << ")";
}