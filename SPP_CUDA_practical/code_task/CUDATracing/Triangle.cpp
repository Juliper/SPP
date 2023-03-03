#include "Triangle.h"

#include "Constants.h"
#include "DirectionPoint.h"
#include "Log.h"
#include "dualray.h"

#include <math.h>

Triangle::Triangle(Point p0, Point p1, Point p2) {
  this->p0 = p0;
  this->p1 = p1;
  this->p2 = p2;
}

bool Triangle::testIntersection(Ray* r, Plane* rpl, Point& nsp, bool find_point) {
  // Drei Kanten berechnen
  DualRay edge1 = DualRay(&(this->p0), &(this->p1));
  DualRay edge2 = DualRay(&(this->p1), &(this->p2));
  DualRay edge3 = DualRay(&(this->p2), &(this->p0));

  float indi1 = r->ray7 * edge1.ray18 + r->ray8 * edge1.ray20 + r->ray9 * edge1.ray22 + r->ray11 * edge1.ray23 +
                r->ray12 * edge1.ray25 + r->ray14 * edge1.ray26;

  float indi2 = r->ray7 * edge2.ray18 + r->ray8 * edge2.ray20 + r->ray9 * edge2.ray22 + r->ray11 * edge2.ray23 +
                r->ray12 * edge2.ray25 + r->ray14 * edge2.ray26;

  float indi3 = r->ray7 * edge3.ray18 + r->ray8 * edge3.ray20 + r->ray9 * edge3.ray22 + r->ray11 * edge3.ray23 +
                r->ray12 * edge3.ray25 + r->ray14 * edge3.ray26;

  int side = -1;

  if (indi1 > 0 && indi2 > 0 && indi3 > 0) {
    side = 1;
  }

  else if (indi1 < 0 && indi2 < 0 && indi3 < 0) {
    side = -1;
  }

  else if (indi1 >= 0 && indi2 >= 0 && indi3 >= 0) {
    side = -1;
  }

  if (side == -1)
    return false;

  // Ansonsten nach Schnitt suchen und dessen Eigenschaften prüfen
  Plane pl = Plane(&(this->p0), &(this->p1), &(this->p2));
  Point sp = Point(&pl, r);
  // FlatPoint fp = FlatPoint(&pl, r);

  float dis2 = ((*rpl) % (&sp));
  // float dis = (*rpl)%(fp);

  // Schranke für den Abstand
  if (dis2 < 2e-3f) {
    return false;
  }

  // falls nur von Interesse ob Schnitt oder nicht kann hier rausgesprungen werden
  if (find_point == false)
    return true;

  // Weitere Details des Punktes berechnen
  Direction u = Direction(p1.point2, p1.point3, p1.point4) - Direction(p0.point2, p0.point3, p0.point4);
  Direction v = Direction(p2.point2, p2.point3, p2.point4) - Direction(p0.point2, p0.point3, p0.point4);

  float uu = u % u;
  float uv = u % v;
  float vv = v % v;

  float D = uv * uv - uu * vv;

  Direction w = Direction(sp.point2, sp.point3, sp.point4) - Direction(p0.point2, p0.point3, p0.point4);
  float wu    = w % u;
  float wv    = w % v;

  float s = (uv * wv - vv * wu) / D;
  float t = (uv * wu - uu * wv) / D;

  // float st = 1-s-t;
  // temp = (st * SCALE);
  sp.n = p0.n * (1 - s - t) + p1.n * s + p2.n * t;
  nsp  = sp;
  // sp.n.normalize();

  return true;
}

bool Triangle::testIntersection(Ray3D* ray) {
  Point orig    = ray->ori;
  Direction dir = ray->dir;

  Direction edge1 = Direction(&p1, &p0);
  Direction edge2 = Direction(&p2, &p0);

  Direction pvec = dir * edge2;
  float det      = edge1 % pvec;

  if (det > -2e-6f && det < 2e-6f)
    return false;

  float inv_det = 1.0f / det;

  Direction tvec = Direction(&orig, &p0);

  float u = (tvec % pvec) * inv_det;
  if (u < 2e-6f || u > 1.00f + 2e-6f)
    return 0;

  Direction qvec = tvec * edge1;

  float v = (dir % qvec) * inv_det;
  if (v < 2e-6f || u + v > 1.0f + 2e-6f)
    return false;

  float t = (edge2 % qvec) * inv_det;
  if (t < 2e-1f)
    return false;

  this->u_t = u;
  this->v_t = v;
  return true;
}

Point Triangle::findIntersectionPoint() {
  float u = this->u_t;
  float v = this->v_t;

  Point ip = Point((1 - u - v) * p0.point2 + u * p1.point2 + v * p2.point2,
                   (1 - u - v) * p0.point3 + u * p1.point3 + v * p2.point3,
                   (1 - u - v) * p0.point4 + u * p1.point4 + v * p2.point4);

  ip.n = Direction((1 - u - v) * p0.n.direction2 + u * p1.n.direction2 + v * p2.n.direction2,
                   (1 - u - v) * p0.n.direction3 + u * p1.n.direction3 + v * p2.n.direction3,
                   (1 - u - v) * p0.n.direction4 + u * p1.n.direction4 + v * p2.n.direction4);
  return ip;
}
