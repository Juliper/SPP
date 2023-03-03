#ifndef DIRECTIONPOINT_H
#define DIRECTIONPOINT_H

#include "Ray.h"
#include "Flatpoint.h"
#include "Plane.h"

class Point;
class Direction;
class Ray3D;
class Ray;
class Camera;
class FlatPoint;

class Direction {
public:
	float direction2;
	float direction3;
	float direction4;

	Direction();
	Direction(float x, float y, float z);
	Direction(Point* p0, Point* p1);
	Direction(Ray* r);
	Direction(float scale, Direction* d);
	Direction(Camera* cam, int x, int y);

	void normalize();
	void invertDirection();
	void printDirection();
	float abs();

	// überladene Operatoren
	bool operator==(Direction d) const;
	Direction operator+ (Direction d) const;
	Direction operator- (Direction d) const;
	float operator%(Direction d) const;
	Direction operator*(Direction d) const;
	Direction operator*(float v) const;
	Direction operator/(float v) const;
};

class Point {
public:
	float point2;
	float point3;
	float point4;
	float point5;
	float point6;
	float u;
	float v;
	Direction n;

	Point ();
	Point (Point* p);
	Point (float x, float y, float z);
	Point (Plane* pl, Ray* r);	
	Point (FlatPoint* fp);

	Point (Point* p, bool d);

	float operator%(Point p) const;
	bool operator==(Point p) const;
	Direction operator- (Point p) const;

	void normalize();
};

class Ray3D
{
public:
	Point ori;
	Direction dir;

	Ray3D();
	Ray3D(Point ori, Direction dir);
	Ray3D(Point p0, Point p1);
};

#endif