#include "MeshSphere.h"

Mesh::Mesh() {
  this->facecount   = 0;
  this->vtxcount    = 0;
  this->objectcount = 0;
}

void Mesh::setFaceCount(int i) {
  faces.resize(i);
  this->facecount = i;
}

void Mesh::setVtxCount(int c) {
  points.resize(c);
  this->vtxcount = c;
}

void Mesh::addFace(int idx, Face f) {
  if ((idx < 0) || (idx >= (int)faces.size()))
    throw std::string("mesh::set_vertex(): index out of range");
  faces[idx] = f;
}

void Mesh::addFace(Face f) {
  this->faces.push_back(f);
  this->facecount++;
}

void Mesh::addPoint(int idx, Point p) {
  if ((idx < 0) || (idx >= (int)points.size()))
    throw std::string("mesh::set_face(): index out of range");
  points[idx] = p;
}

void Mesh::addPoint(Point p) {
  this->points.push_back(p);
  this->vtxcount++;
}

void Mesh::configureMesh() {
  for (int i = 0; i < this->facecount; i++) {
    Triangle t = Triangle(points[faces[i].p0], points[faces[i].p1], points[faces[i].p2]);
    this->addObject(t);
  }

  this->bound = Sphere(this);
}

void Mesh::addObject(Triangle t) {
  objects.push_back(t);
  this->objectcount++;
}

int Mesh::getObjectCount() {
  return this->objectcount;
}

void Mesh::interpolateNormals() {
  // Flächennormalen berechnen
  for (int i = 0; i < this->facecount; i++) {
    // Kanten des Dreiecks
    Direction edge1 = points[faces[i].p0] - points[faces[i].p1];
    Direction edge2 = points[faces[i].p1] - points[faces[i].p2];

    // Flächennormale durch Kreuzprodukt
    Direction fnorm = edge1 * edge2;
    fnorm.normalize();
    this->facenormals.push_back(fnorm);
  }

  // Für jeden Punkt Normale aus den Flächennormalen interpolieren an denen der Punkt beteiligt ist
  for (int i = 0; i < this->vtxcount; i++) {
    int cnt       = 0;
    Direction sum = Direction(0.0f, 0.0f, 0.0f);

    for (int j = 0; j < this->facecount; j++) {
      if (faces[j].p0 == i || faces[j].p1 == i || faces[j].p2 == i) {
        cnt++;
        sum = sum + this->facenormals[j];
      }
    }

    points[i].n = (sum / (float)cnt);
    points[i].n.normalize();
  }
}
