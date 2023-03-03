#include "Pointpair.h"

#include <math.h>
#include <stdio.h>

double DistanceBetweenPoints(Point* p0, Point* p1) {
  // return sqrt(abs(p0->point2 * p1->point2 + p0->point3 * p1->point3 + p0->point4 * p1->point4));
  return sqrt((p0->point2 - p1->point2) * (p0->point2 - p1->point2) +
              (p0->point3 - p1->point3) * (p0->point3 - p1->point3) +
              (p0->point4 - p1->point4) * (p0->point4 - p1->point4));
}

PointPair::PointPair(Ray* r, Sphere* s) {
  this->PP7  = -r->ray14;
  this->PP8  = r->ray12;
  this->PP9  = -s->s3 * r->ray14 - r->ray11 * s->s5 + s->s4 * r->ray12;
  this->PP10 = r->ray11;
  this->PP11 = -r->ray9;
  this->PP12 = -s->s4 * r->ray9 + r->ray8 * s->s5 + s->s2 * r->ray14;
  this->PP13 = -r->ray8;
  this->PP14 = -r->ray7 * s->s5 - s->s2 * r->ray12 + s->s3 * r->ray9;
  this->PP15 = r->ray7;
  this->PP16 = -s->s3 * r->ray8 + s->s4 * r->ray7 + s->s2 * r->ray11;

  this->inprod = 2 * PP12 * PP13 - PP11 * PP11 - PP7 * PP7 + 2 * PP14 * PP15 + 2 * PP9 * PP10 + PP16 * PP16 - PP8 * PP8;
}

bool PointPair::isPPReal() {
  if (inprod >= 2e-5)
    return true;
  else
    return false;
}

/*
Momentan wird jetzt nur pminus ausgeben. Es scheint, das automatisch immer der Punkt mit geringstem Abstand gewählt wird
*/
Point* PointPair::extractPoint() {
  float w       = sqrt(this->inprod);
  Point* pminus = new Point(0, 0, 0);

  pminus->point2 = -PP8 * PP15 - PP7 * PP13 + PP16 * PP10 + w * PP10;
  pminus->point3 = -PP11 * PP15 + PP7 * PP10 + w * PP13 + PP16 * PP13;
  pminus->point4 = PP8 * PP10 + PP16 * PP15 + PP11 * PP13 + w * PP15;
  pminus->point5 = w * PP16 + PP13 * PP12 + PP15 * PP14 + PP10 * PP9 + PP16 * PP16;
  pminus->point6 = PP15 * PP15 + PP10 * PP10 + PP13 * PP13;
  return pminus;
}
