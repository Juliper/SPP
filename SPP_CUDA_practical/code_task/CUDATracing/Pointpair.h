#ifndef POINTPAIR_H
#define POINTPAIR_H

#include "Ray.h"
#include "MeshSphere.h"
#include "DirectionPoint.h"

class Sphere;

class PointPair {
public:
	float PP7;
	float PP8;
	float PP9;
	float PP10;
	float PP11;
	float PP12;
	float PP13;
	float PP14;
	float PP15;
	float PP16;
	float inprod;

	PointPair(Ray* r, Sphere* s);
	bool isPPReal();
	bool isPPReal(int recu);
	Point * extractPoint();
};

#endif