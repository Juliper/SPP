#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "DirectionPoint.h"
#include "Plane.h"

class Triangle {
public:
	Point p0;
	Point p1;
	Point p2;

	float u_t;
	float v_t;

	Triangle (Point p0, Point p1, Point p2);
	Triangle ()	{}
	bool testIntersection(Ray* r, Plane* rpl, Point& nsp, bool find_point);
	bool testIntersection(Ray3D* ray);
	Point findIntersectionPoint();
};

#endif