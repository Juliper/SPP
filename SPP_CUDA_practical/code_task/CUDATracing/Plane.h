#ifndef PLANE_H
#define PLANE_H

#include "DirectionPoint.h"
#include "Ray.h"
#include "Flatpoint.h"

class Direction;
class Point;
class Ray;
class FlatPoint;

class Plane {
public:
	float pl2;
	float pl3;
	float pl4;
	float pl5;

	Plane(Direction* n, float dinstance);
	Plane(Ray* r);
	Plane(float n1, float n2, float n3, float d);
	Plane(Point* p0, Point* p1, Point* p2);

	
	float operator%(Point* p) const;
	float operator%(FlatPoint fp) const;
};
#endif