#include "Ray3D.h"

Ray3D::Ray3D() {
}

Ray3D::Ray3D(Point ori, Direction dir) {
  this->ori = ori;
  this->dir = dir;
}

Ray3D::Ray3D(Point p0, Point p1) {
  this->ori      = p0;
  Direction temp = Direction(&p1, &p0);
  temp.normalize();
  this->dir = temp;
}