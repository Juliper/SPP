#include "Plane.h"

#include <stdio.h>

Plane::Plane(Direction* n, float distance) {
  this->pl2 = n->direction2;
  this->pl3 = n->direction3;
  this->pl4 = n->direction4;
  this->pl5 = distance;
}

Plane::Plane(Ray* r) {
  Direction* n = r->getDirection();
  // n->normalize();

  Point* o = r->origin;

  this->pl2 = n->direction2;
  this->pl3 = n->direction3;
  this->pl4 = n->direction4;
  this->pl5 = (o->point2 * n->direction2 + o->point3 * n->direction3 + o->point4 * n->direction4) / o->point6;
  delete n;
}

Plane::Plane(Point* p0, Point* p1, Point* p2) {
  this->pl2 = -p2->point3 * p0->point4 + p2->point4 * p0->point3 + p0->point4 * p1->point3 + p2->point3 * p1->point4 -
              p2->point4 * p1->point3 - p0->point3 * p1->point4;

  this->pl3 = -p2->point4 * p0->point2 - p0->point4 * p1->point2 + p0->point2 * p1->point4 + p2->point2 * p0->point4 +
              p2->point4 * p1->point2 - p2->point2 * p1->point4;

  this->pl4 = p0->point3 * p1->point2 * p2->point6 - p0->point2 * p1->point3 * p2->point6 -
              p2->point3 * p1->point2 * p0->point6 + p2->point3 * p0->point2 * p1->point6 -
              p2->point2 * p0->point3 * p1->point6 + p2->point2 * p1->point3 * p0->point6;

  this->pl5 = p2->point3 * p0->point2 * p1->point4 + p2->point4 * p0->point3 * p1->point2 +
              p2->point2 * p0->point4 * p1->point3 - p2->point4 * p0->point2 * p1->point3 -
              p2->point2 * p0->point3 * p1->point4 - p2->point3 * p0->point4 * p1->point2;
}

Plane::Plane(float n1, float n2, float n3, float d) {
  this->pl2 = n1;
  this->pl3 = n2;
  this->pl4 = n3;
  this->pl5 = d;
}

float Plane::operator%(Point* p) const {
  return (pl2 * p->point2 + pl3 * p->point3 + pl4 * p->point4) - pl5 * p->point6;
}

float Plane::operator%(FlatPoint fp) const {
  return -pl2 * fp.fp9 - pl3 * fp.fp12 - pl4 * fp.fp14 - pl5 * fp.fp16;
}
