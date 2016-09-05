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
	glClearColor(0.9,0.9,0.9,1.0);
}

void RenderWidget::draw()
{
	const float nbSteps = 200.0;

	glBegin(GL_QUAD_STRIP);
	for (int i = 0; i<nbSteps; ++i)
	{
		const float ratio = i / nbSteps;
		const float angle = 21.0*ratio;
		const float c = cos(angle);
		const float s = sin(angle);
		const float r1 = 1.0 - 0.8f*ratio;
		const float r2 = 0.8f - 0.8f*ratio;
		const float alt = ratio - 0.5f;
		const float nor = 0.5f;
		const float up = sqrt(1.0 - nor*nor);
		glColor3f(1.0 - ratio, 0.2f, ratio);
		glNormal3f(nor*c, up, nor*s);
		glVertex3f(r1*c, alt, r1*s);
		glVertex3f(r2*c, alt + 0.05f, r2*s);
	}
	glEnd();
}