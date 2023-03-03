#include "Ray.h"

#include <math.h>
#include <stdio.h>

Ray::Ray(Camera* cam, int x, int y) {
  // Zuerst Direction des Strahles erstellen
  Direction* dir = new Direction(cam, x, y);

  // Hier müsste pos6 eigentlich jeweils 1 sein! Daher die Multiplikation in 7,8 und 11 entfernt
  this->ray7  = -dir->direction4;
  this->ray8  = dir->direction3;
  this->ray9  = cam->pos->point4 * dir->direction3 - cam->pos->point3 * dir->direction4;
  this->ray11 = -dir->direction2;
  this->ray12 = -cam->pos->point4 * dir->direction2 + cam->pos->point2 * dir->direction4;
  this->ray14 = cam->pos->point3 * dir->direction2 - cam->pos->point2 * dir->direction3;

  this->origin = cam->pos;
  delete dir;
}

Ray::Ray(Point* p0, Point* p1) {
  // alle p6 == 1 angenommen, da Punkte normalisiert
  this->ray7   = p0->point4 - p1->point4;
  this->ray8   = p1->point3 - p0->point3;
  this->ray9   = p0->point4 * p1->point3 - p0->point3 * p1->point4;
  this->ray11  = -p1->point2 + p0->point2;
  this->ray12  = -p0->point4 * p1->point2 + p0->point2 * p1->point4;
  this->ray14  = p0->point3 * p1->point2 - p0->point2 * p1->point3;
  this->origin = p0;
}

Ray::Ray(Point* p0, Direction* dir) {
  this->ray7  = -dir->direction4;
  this->ray8  = dir->direction3;
  this->ray9  = p0->point4 * dir->direction3 - p0->point3 * dir->direction4;
  this->ray11 = -dir->direction2;
  this->ray12 = p0->point2 * dir->direction4 - p0->point4 * dir->direction2;
  this->ray14 = p0->point3 * dir->direction2 - p0->point2 * dir->direction3;

  this->origin = p0;
}

Ray::Ray() {
  this->ray7  = 0;
  this->ray8  = 0;
  this->ray9  = 0;
  this->ray11 = 0;
  this->ray12 = 0;
  this->ray14 = 0;
}

Direction* Ray::getDirection() {
  return new Direction(-this->ray11, this->ray8, -this->ray7);
}

float Ray::abs() {
  return (float)sqrt(ray7 * ray7 + ray8 * ray8 + ray11 * ray11);
}

void Ray::normalize() {
  float absval = this->abs();

  this->ray7  = ray7 / absval;
  this->ray8  = ray8 / absval;
  this->ray9  = ray9 / absval;
  this->ray11 = ray11 / absval;
  this->ray12 = ray12 / absval;
  this->ray14 = ray14 / absval;
}