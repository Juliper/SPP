#ifndef MATERIAL_H
#define MATERIAL_H
#include "Color.h"
#include "Constants.h"

class Material {
public:
	Color ambient;
	Color diffuse;
	Color specular;
	float shininess;

	Material();
	Material(Color a, Color d, Color s, float ex);
};

#endif