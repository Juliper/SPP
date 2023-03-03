#ifndef FLATPOINT_H
#define FLATPOINT_H

#include "DirectionPoint.h"
#include "Plane.h"
#include "Ray.h"

class Ray;
class Plane;

class FlatPoint {
public:
	float fp9;
	float fp12;
	float fp14;
	float fp16;

	FlatPoint(Plane* p, Ray* r);
};
#endif