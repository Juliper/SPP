#include "Transformation.h"

#include "math.h"

Transformation::Transformation() {
  for (int i = 0; i < 15; i++) {
    this->m[i] = 0.0;
  }

  this->m[0]  = 1.0f;
  this->m[5]  = 1.0f;
  this->m[10] = 1.0f;
  this->m[15] = 1.0f;
}

void Transformation::setRotation(float angle, float x, float y, float z) {
  Direction axis(x, y, z);
  axis.normalize();

  float b = angle;  // / 3.141592653589793f;

  float c  = cosf(b);
  float ac = 1.00f - c;
  float s  = sinf(b);

  m[0] = axis.direction2 * axis.direction2 * ac + c;
  m[1] = axis.direction2 * axis.direction3 * ac + axis.direction4 * s;
  m[2] = axis.direction2 * axis.direction4 * ac - axis.direction3 * s;

  m[4] = axis.direction3 * axis.direction2 * ac - axis.direction4 * s;
  m[5] = axis.direction3 * axis.direction3 * ac + c;
  m[6] = axis.direction3 * axis.direction4 * ac + axis.direction2 * s;

  m[8]  = axis.direction4 * axis.direction2 * ac + axis.direction3 * s;
  m[9]  = axis.direction4 * axis.direction3 * ac - axis.direction2 * s;
  m[10] = axis.direction4 * axis.direction4 * ac + c;
}

void Transformation::setTranslation(float x, float y, float z) {
  this->m[12] = x;
  this->m[13] = y;
  this->m[14] = z;
  this->m[15] = 1.0;
}

Point Transformation::transformPoint(Point p) {
  float xn = p.point2 * m[0] + p.point3 * m[4] + p.point4 * m[8] + m[12];
  float yn = p.point2 * m[1] + p.point3 * m[5] + p.point4 * m[9] + m[13];
  float zn = p.point2 * m[2] + p.point3 * m[6] + p.point4 * m[10] + m[14];
  float wn = p.point2 * m[3] + p.point3 * m[7] + p.point4 * m[11] + m[15];

  float xnn = p.n.direction2 * m[0] + p.n.direction3 * m[4] + p.n.direction4 * m[8];   //+ m[12];
  float ynn = p.n.direction2 * m[1] + p.n.direction3 * m[5] + p.n.direction4 * m[9];   //+ m[13];
  float znn = p.n.direction2 * m[2] + p.n.direction3 * m[6] + p.n.direction4 * m[10];  // + m[14];

  Point pnew     = Point(xn / wn, yn / wn, zn / wn);
  Direction nnew = Direction(xnn, ynn, znn);
  pnew.n         = nnew;
  nnew.normalize();
  pnew.u = p.u;
  pnew.v = p.v;

  return pnew;
}