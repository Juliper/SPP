#include "MeshSphere.h"
#include "Plane.h"

#include <algorithm>
#include <float.h>
#include <math.h>
#include <stdio.h>

Sphere::Sphere() {
}

Sphere::Sphere(Point p, float radius) {
  this->radius = radius;
  this->MP     = p;

  this->s2 = p.point2;
  this->s3 = p.point3;
  this->s4 = p.point4;
  this->s5 = p.point5 - 0.5f * radius * radius;
  this->s6 = p.point6;
}

Sphere::Sphere(float s2, float s3, float s4, float s5) {
  this->s2 = s2;
  this->s3 = s3;
  this->s4 = s4;
  this->s5 = s5;
  this->s6 = 1.0f;
}

Sphere::Sphere(Mesh* m) {
  Plane dP1 = Plane(-1, 0, 0, 0);
  Plane dP2 = Plane(0, -1, 0, 0);
  Plane dP3 = Plane(0, 0, -1, 0);

  float minDis[3] = {1e10, 1e10, 1e10};
  float maxDis[3] = {-1e10, -1e10, -1e10};

  // Abstand für jede Koordinate jedes Punktes zu den drei Ebenen berechnen
  for (int i = 0; i < m->vtxcount; i++) {
    Point p = m->points[i];

    // Distanz zum ersten Punkt p0 des Dreieicks
    minDis[0] = min(p.point2, minDis[0]);
    maxDis[0] = max(p.point2, maxDis[0]);

    // Distanz zum zweiten Punkt p1 des Dreieicks
    minDis[1] = min(p.point3, minDis[1]);
    maxDis[1] = max(p.point3, maxDis[1]);

    // Distanz zum zweiten Punkt p1 des Dreieicks
    minDis[2] = min(p.point4, minDis[2]);
    maxDis[2] = max(p.point4, maxDis[2]);
  }

  // Mittelpunkt der Bounding Sphere
  Point MP = Point(0.5f * (minDis[0] + maxDis[0]), 0.5f * (minDis[1] + maxDis[1]), 0.5f * (minDis[2] + maxDis[2]));

  // Radius der Bounding Sphere
  float r = 0.0f;
  for (int i = 0; i < m->vtxcount; i++) {
    float lan = MP % (m->points[i]);
    r         = min(lan, r);
  }

  this->radius = sqrt(abs(2.0f * r - 2.0f * MP.point5));
  this->MP     = &MP;

  this->s2 = MP.point2;
  this->s3 = MP.point3;
  this->s4 = MP.point4;
  this->s5 = r;
  this->s6 = MP.point6;
}

Sphere::Sphere(Mesh2* m, std::vector<Point> points) {
  Plane dP1 = Plane(-1, 0, 0, 0);
  Plane dP2 = Plane(0, -1, 0, 0);
  Plane dP3 = Plane(0, 0, -1, 0);

  float minDis[3] = {1e10, 1e10, 1e10};
  float maxDis[3] = {-1e10, -1e10, -1e10};

  // Abstand für jede Koordinate jedes Punktes zu den drei Ebenen berechnen
  for (int i = 0; i < m->vtxcount; i++) {
    Point p = points[m->points[i]];

    // Distanz zum ersten Punkt p0 des Dreieicks
    minDis[0] = min(p.point2, minDis[0]);
    maxDis[0] = max(p.point2, maxDis[0]);

    // Distanz zum zweiten Punkt p1 des Dreieicks
    minDis[1] = min(p.point3, minDis[1]);
    maxDis[1] = max(p.point3, maxDis[1]);

    // Distanz zum zweiten Punkt p1 des Dreieicks
    minDis[2] = min(p.point4, minDis[2]);
    maxDis[2] = max(p.point4, maxDis[2]);
  }

  // Mittelpunkt der Bounding Sphere
  Point MP = Point(0.5f * (minDis[0] + maxDis[0]), 0.5f * (minDis[1] + maxDis[1]), 0.5f * (minDis[2] + maxDis[2]));

  // Radius der Bounding Sphere
  float r = 0.0f;
  for (int i = 0; i < m->vtxcount; i++) {
    float lan = MP % (points[m->points[i]]);
    r         = min(lan, r);
  }

  this->radius = sqrt(abs(2.0f * r - 2.0f * MP.point5));
  this->MP     = &MP;

  this->s2 = MP.point2;
  this->s3 = MP.point3;
  this->s4 = MP.point4;
  this->s5 = r;
  this->s6 = MP.point6;
}

Point* Sphere::giveMP() {
  return &(this->MP);
}

bool Sphere::doesRayIntersect(Ray* r) {
  PointPair PP = PointPair(r, this);
  return PP.isPPReal();
}

bool Sphere::doesRayIntersect(Ray3D* r) {
  // r->dir.normalize();
  Direction oc = Direction(&(r->ori), &(MP));

  double b = r->dir % (oc);
  double c = (oc % oc) - (this->radius * this->radius);

  if (b * b - c < 0.00001)
    return false;
  else
    return true;
}
bool Sphere::doesRayIntersect2(Ray* r) {
  float inprod = 2.0f * (-r->ray11 * s3 * r->ray14 - r->ray7 * r->ray7 * s5 + s3 * r->ray9 * r->ray7 +
                         r->ray8 * s4 * r->ray9 - r->ray8 * r->ray8 * s5 + s4 * r->ray7 * s2 * r->ray11 -
                         s3 * r->ray8 * s4 * r->ray7 - s2 * r->ray12 * r->ray7 + s4 * r->ray12 * r->ray11 -
                         s2 * r->ray14 * r->ray8 - s3 * r->ray8 * s2 * r->ray11 - r->ray11 * r->ray11 * s5)

                 + s2 * s2 * r->ray11 * r->ray11 - r->ray9 * r->ray9 + s3 * s3 * r->ray8 * r->ray8 -
                 r->ray12 * r->ray12 - r->ray14 * r->ray14 + s4 * s4 * r->ray7 * r->ray7;
  if (inprod >= 0.0f)
    return true;

  return false;
}

float Sphere::distOfPoints(Point* p1, Point* p2) {
  return sqrt((p1->point2 - p2->point2) * (p1->point2 - p2->point2) +
              (p1->point3 - p2->point3) * (p1->point3 - p2->point3) +
              (p1->point4 - p2->point4) * (p1->point4 - p2->point4));
}