#include "dualray.h"

DualRay::DualRay(Point* p0, Point* p1) {
  this->ray18 = -p0->point3 * p1->point2 + p0->point2 * p1->point3;  // e1^e2^einf;
  this->ray20 = -p0->point4 * p1->point2 + p0->point2 * p1->point4;  // e1^e3^einf;
  this->ray22 = p1->point2 * p0->point6 - p0->point2 * p1->point6;   // e1^einf^e0
  this->ray23 = p0->point3 * p1->point4 - p0->point4 * p1->point3;   // e2^e3^einf
  this->ray25 = p1->point3 * p0->point6 - p0->point3 * p1->point6;   // e2^einf^e0
  this->ray26 = p1->point4 * p0->point6 - p0->point4 * p1->point6;   // e3^einf^e0
}