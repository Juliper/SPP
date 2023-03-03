#ifndef LACPUTRACERGLUT_H
#define LACPUTRACERGLUT_H

#include "virtualTracer_glut.h"
#include "Constants.h"
#include "Ray.h"
#include "DirectionPoint.h"
#include "Camera.h"
#include "Plane.h"
#include "Color.h"
#include "Light.h"

class LA_CPU_Tracer : public virtualTracer
{
public:
	LA_CPU_Tracer(int width, int height);

	Direction reflect3DLA(Direction* l, Direction* n);
	float calculateDistanceBetweenPoints3D(Point* p0, Point* p1);
	Color traceRay(Ray3D* ray ,Camera* cam, Color col, int recu);
	void rayTraceScene(void);

	int initialize();
	int run();
	void renderScene();	
	void renderImage();

private:
	float la_cpu_stat;
//	int testcnt;
};
#endif