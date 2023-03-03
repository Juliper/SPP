#ifndef LIGHT_H
#define LIGHT_H

#include "DirectionPoint.h"
#include "Color.h"

class Light {
public:
	Color* ambient;
	Color* diffuse;
	Color* specular;

	Point pos;

	Light(Point pos, Color* a, Color* d, Color* s);
	Light();
	
};

#endif