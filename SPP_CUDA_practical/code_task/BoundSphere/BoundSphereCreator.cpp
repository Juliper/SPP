#include "MeshSphere.h"
#include "Miniball.h"
#include "PointN.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

BoundSphereCreator::BoundSphereCreator() {
}

float* BoundSphereCreator::createBound(Mesh2* m, std::vector<Point> vtx) {
  using std::cout;
  using std::endl;
  const int d = 3;

  PointN p(d);
  Miniball mb(d);

  // Copy Points To Local Format
  for (int i = 0; i < m->vtxcount; ++i) {
    p[0] = vtx[m->points[i]].point2;
    p[1] = vtx[m->points[i]].point3;
    p[2] = vtx[m->points[i]].point4;

    mb.check_in(p);
  }

  // construct ball
  mb.build();

  float* res = (float*)malloc(4 * sizeof(float));
  PointN pc  = mb.center();
  res[0]     = (float)pc[0];
  res[1]     = (float)pc[1];
  res[2]     = (float)pc[2];
  res[3]     = (float)mb.squared_radius();

  return res;
}
