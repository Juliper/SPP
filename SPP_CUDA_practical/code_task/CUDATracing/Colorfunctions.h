#ifndef COLORFUNCTIONS_H
#define COLORFUNCTIONS_H

#include "Material.h"
#include "Light.h"
#include "DirectionPoint.h"

float angleBetweenDirections(Direction* v1, Direction* v2);
Color calculateAmbient(Light* l, Material* mat);
Color calculateDiffuse(Light* l, Material* mat, float NDotL);
Color calculateSpecular(Light* l, Material* mat, float RDotE);

#endif