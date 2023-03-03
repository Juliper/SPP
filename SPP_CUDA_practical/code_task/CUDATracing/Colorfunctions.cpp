#include "Colorfunctions.h"

#include <math.h>
#include <stdio.h>

float angleBetweenDirections(Direction* v1, Direction* v2) {
  v1->normalize();
  v2->normalize();
  float r1DOTr2 = (*v1) % (*v2);

  if (r1DOTr2 < 0)
    r1DOTr2 = 0.0;

  return (float)fabs(r1DOTr2);
}

Color calculateAmbient(Light* l, Material* mat) {
  return Color(*(l->ambient) * mat->ambient);
}

Color calculateDiffuse(Light* l, Material* mat, float NDotL) {
  Color temp = *(l->diffuse) * mat->diffuse;
  Color res  = temp * NDotL;

  return res;
}

Color calculateSpecular(Light* l, Material* mat, float RDotE) {
  Color temp = *(l->specular) * mat->specular;
  return (temp * (float)pow(RDotE, mat->shininess));
}