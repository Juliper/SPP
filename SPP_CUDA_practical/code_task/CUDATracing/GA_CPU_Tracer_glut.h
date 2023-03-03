#ifndef GACPUTRACERGLUT_H
#define GACPUTRACERGLUT_H

#include "virtualTracer_glut.h"
#include "Constants.h"
#include "Ray.h"
#include "DirectionPoint.h"
#include "Camera.h"
#include "Plane.h"
#include "Color.h"
#include "Light.h"

class GA_CPU_Tracer : public virtualTracer
{

public:

	GA_CPU_Tracer(int width, int height);

	Direction reflect3D(Direction* l, Direction* n);
	Ray reflect5D(Ray* l, Direction* n);
	float calculateDistanceBetweenPoints(Point* p0, Point* p1);
	Color traceRay(Ray* ray, Plane* rpl, Camera * cam, Color col, int recursion);
	void rayTraceSceneNew(void);

	int initialize();
	int run();
	void renderScene();	
	void renderImage();

private:
	float ga_cpu_stat;
//	int testcnt;
};
#endif