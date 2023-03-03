#include "Material.h"

Material::Material(Color a, Color d, Color s, float ex) {
  this->ambient   = a;
  this->diffuse   = d;
  this->specular  = s;
  this->shininess = ex;
}

Material::Material() {
}
