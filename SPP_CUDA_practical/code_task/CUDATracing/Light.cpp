#include "Light.h"

Light::Light(Point pos, Color* a, Color* d, Color* s) {
  this->pos      = pos;
  this->ambient  = a;
  this->diffuse  = d;
  this->specular = s;
}

Light::Light() {
}
