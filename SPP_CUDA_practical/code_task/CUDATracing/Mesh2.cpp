#include "MeshSphere.h"

Mesh2::Mesh2() {
  this->facecount = 0;
  this->faces.resize(0);
  this->vtxcount = 0;
  this->points.resize(0);
  this->objectcount = 0;
  this->objects.resize(0);
}

void Mesh2::addFace(Face f) {
  this->faces.push_back(f);
  this->facecount++;
}

void Mesh2::addObject(Triangle t) {
  this->objects.push_back(t);
  this->objectcount++;
}

void Mesh2::addPoints(Face* f) {
  if (pointIsNew(f->p0)) {
    this->points.push_back(f->p0);
    this->vtxcount++;
  }

  if (pointIsNew(f->p1)) {
    this->points.push_back(f->p1);
    this->vtxcount++;
  }

  if (pointIsNew(f->p2)) {
    this->points.push_back(f->p2);
    this->vtxcount++;
  }
}

bool Mesh2::pointIsNew(int p) {
  bool found = true;
  for (int i = 0; i < this->vtxcount; i++) {
    if (this->points[i] == p)
      found = false;
  }
  return found;
}

void Mesh2::configureMesh(std::vector<Point> points) {
  this->objectcount = 0;
  this->objects.resize(0);
  for (int i = 0; i < this->facecount; i++) {
    Triangle t = Triangle(points[faces[i].p0], points[faces[i].p1], points[faces[i].p2]);
    this->addObject(t);
  }

  BoundSphereCreator* bsc = new BoundSphereCreator();
  float* res              = bsc->createBound(this, points);
  float r                 = sqrt(res[3]);
  Point c                 = Point(res[0], res[1], res[2]);

  this->bound = Sphere(c, r);
  // this->bound = Sphere(this, points);
}

int Mesh2::getObjectCount() {
  return this->objectcount;
}
