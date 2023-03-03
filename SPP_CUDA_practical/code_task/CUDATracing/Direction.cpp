#include "DirectionPoint.h"

#include <math.h>
#include <stdio.h>

Direction::Direction() {
  this->direction2 = 0;
  this->direction3 = 0;
  this->direction4 = 0;
}

Direction::Direction(float x, float y, float z) {
  this->direction2 = x;
  this->direction3 = y;
  this->direction4 = z;
}

Direction::Direction(Camera* cam, int x, int y) {
  this->direction2 = x * cam->pixelWidth - 0.5f * cam->cameraFOVwidth;
  this->direction3 = y * cam->pixelHeight - 0.5f * cam->cameraFOVheight;
  this->direction4 = 1;
}

Direction::Direction(Point* p0, Point* p1) {
  this->direction2 = p0->point2 - p1->point2;
  this->direction3 = p0->point3 - p1->point3;
  this->direction4 = p0->point4 - p1->point4;
}

Direction::Direction(Ray* r) {
  this->direction2 = -r->ray11;
  this->direction3 = r->ray8;
  this->direction4 = -r->ray7;
}

Direction::Direction(float scale, Direction* d) {
  this->direction2 = scale * d->direction2;
  this->direction3 = scale * d->direction3;
  this->direction4 = scale * d->direction4;
}

void Direction::normalize() {
  float absval     = abs();
  this->direction2 = direction2 / absval;
  this->direction3 = direction3 / absval;
  this->direction4 = direction4 / absval;
}

void Direction::invertDirection() {
  this->direction2 = -direction2;
  this->direction3 = -direction3;
  this->direction4 = -direction4;
}

float Direction::abs() {
  return (float)sqrt(direction2 * direction2 + direction3 * direction3 + direction4 * direction4);
}

void Direction::printDirection() {
  printf("Direction e1:%f e2:%f e3:%f \n", direction2, direction3, direction4);
}

bool Direction::operator==(Direction d) const {
  return this->direction2 == d.direction2 && direction3 == d.direction3 && direction4 == d.direction4;
}

Direction Direction::operator+(Direction d) const {
  return Direction(direction2 + d.direction2, direction3 + d.direction3, direction4 + d.direction4);
}

Direction Direction::operator-(Direction d) const {
  return Direction(direction2 - d.direction2, direction3 - d.direction3, direction4 - d.direction4);
}

float Direction::operator%(Direction d) const {
  return direction2 * d.direction2 + direction3 * d.direction3 + direction4 * d.direction4;
}

Direction Direction::operator*(Direction d) const {
  return Direction(direction3 * d.direction4 - direction4 * d.direction3,
                   direction4 * d.direction2 - direction2 * d.direction4,
                   direction2 * d.direction3 - direction3 * d.direction2);
}

Direction Direction::operator*(float v) const {
  return Direction(direction2 * v, direction3 * v, direction4 * v);
}

Direction Direction::operator/(float v) const {
  return Direction(direction2 / v, direction3 / v, direction4 / v);
}