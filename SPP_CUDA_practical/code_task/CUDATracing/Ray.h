#ifndef RAY_H
#define RAY_H

#include "Camera.h"
#include "DirectionPoint.h"
//#include "point.h"

class Direction;
class Point;
class Camera;

class Ray {
public:
	float ray7;
	float ray8;
	float ray9;
	float ray11;
	float ray12;
	float ray14;
	Point* origin;

	Ray(Camera* cam, int x, int y);
	Ray(Point* p0, Point* p1);
	Ray(Point* p0, Direction* dir);
	Ray();
	Direction* getDirection();
	float abs();
	void normalize();
};

#endif