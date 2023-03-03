#include "DirectionPoint.h"

#include <math.h>
#include <stdio.h>

Point::Point() {
  this->point2 = 0;
  this->point3 = 0;
  this->point4 = 0;
  this->point5 = 0;
  this->point6 = 0;
}

Point::Point(Point* p) {
  this->point2 = p->point2;
  this->point3 = p->point3;
  this->point4 = p->point4;
  this->point5 = p->point5;
  this->point6 = p->point6;
}

Point::Point(float x, float y, float z) {
  this->point2 = x;
  this->point3 = y;
  this->point4 = z;
  this->point5 = 0.5f * (x * x + y * y + z * z);
  this->point6 = 1.0f;
}

Point::Point(FlatPoint* fp) {
  this->point2 = fp->fp9;
  this->point3 = fp->fp12;
  this->point4 = fp->fp14;
  this->point5 = 0.5f * (fp->fp9 * fp->fp9 + fp->fp12 * fp->fp12 + fp->fp14 * fp->fp14);
  this->point6 = fp->fp16;
}

Point::Point(Plane* pl, Ray* r) {
  FlatPoint* fp = new FlatPoint(pl, r);

  this->point2 = fp->fp9;
  this->point3 = fp->fp12;
  this->point4 = fp->fp14;
  this->point5 = 0.5f * (fp->fp9 * fp->fp9 + fp->fp12 * fp->fp12 + fp->fp14 * fp->fp14);
  this->point6 = fp->fp16;
  delete fp;
}

float Point::operator%(Point p) const {
  return (this->point2 * p.point2 + this->point3 * p.point3 + this->point4 * p.point4) - this->point5 * p.point6 -
         this->point6 * p.point5;
}

bool Point::operator==(Point p) const {
  return (this->point2 == p.point2 && this->point3 == p.point3 && this->point4 == p.point4);
}

Direction Point::operator-(Point p) const {
  return Direction(point2 / point6 - p.point2 / p.point6, point3 / point6 - p.point3 / p.point6,
                   point4 / point6 - p.point4 / p.point6);
}