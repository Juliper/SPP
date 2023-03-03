#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "DirectionPoint.h"

class Transformation
{
private:
	float m[16];

public:
	Transformation();
	void setTranslation(float x, float y, float z);
	void setRotation(float angle, float x, float y, float z);
	void setScaling(float f);
	
	Point transformPoint(Point p);
};

#endif