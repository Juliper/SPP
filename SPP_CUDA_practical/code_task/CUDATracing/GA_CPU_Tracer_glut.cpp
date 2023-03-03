#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <windows.h>

#include "Camera.h"
#include "Color.h"
#include "Colorfunctions.h"
#include "Constants.h"
#include "DirectionPoint.h"
#include "GA_CPU_Tracer_glut.h"
#include "Light.h"
#include "Log.h"
#include "Material.h"
#include "MeshSphere.h"
#include "Plane.h"
#include "Ray.h"
#include "Triangle.h"
#include "timing.h"

#include <omp.h>

GA_CPU_Tracer::GA_CPU_Tracer(int width, int height) : virtualTracer(width, height) {
  this->width       = width;
  this->height      = height;
  this->checkImage  = new float[width * height * 3];
  this->ga_cpu_stat = 0;
}

Direction GA_CPU_Tracer::reflect3D(Direction* l, Direction* n) {
  n->normalize();
  l->normalize();

  Direction dir_ref = Direction(0.1f, 0, 0);

  // 1. Stufe
  double n2Q = n->direction2 * n->direction2;
  double n3Q = n->direction3 * n->direction3;
  double n4Q = n->direction4 * n->direction4;

  // 2. Stufe
  dir_ref.direction2 =
      (float)(2.0f * n->direction2 * n->direction3 * l->direction3 +
              2.0f * n->direction2 * n->direction4 * l->direction4 + l->direction2 * (n2Q - n4Q - n3Q));

  dir_ref.direction3 = (float)(2 * n->direction3 * n->direction2 * l->direction2 +
                               2 * n->direction3 * n->direction4 * l->direction4 + l->direction3 * (n3Q - n2Q - n4Q));

  dir_ref.direction4 = (float)(2 * n->direction4 * n->direction2 * l->direction2 +
                               2 * n->direction4 * n->direction3 * l->direction3 + l->direction4 * (n4Q - n2Q - n3Q));
  return dir_ref;
}

Ray GA_CPU_Tracer::reflect5D(Ray* l, Direction* n) {
  Point p  = l->origin;
  double d = n->direction2 * p.point2 + n->direction3 * p.point3 + n->direction4 * p.point4;

  Plane pl = Plane(n, (float)d);
  Ray r    = Ray();
  r.origin = l->origin;

  // 1. Schritt
  double p2Q = pl.pl2 * pl.pl2;
  double p3Q = pl.pl3 * pl.pl3;
  double p4Q = pl.pl4 * pl.pl4;

  double two23 = 2.0f * pl.pl2 * pl.pl3;
  double two24 = 2.0f * pl.pl2 * pl.pl4;
  double two34 = 2.0f * pl.pl3 * pl.pl4;

  double two25 = 2.0f * pl.pl2 * pl.pl5;
  double two35 = 2.0f * pl.pl3 * pl.pl5;
  double two45 = 2.0f * pl.pl4 * pl.pl5;

  // 2. Schritt
  r.ray7 = (float)(two24 * l->ray11 - two34 * l->ray8 + l->ray7 * (p4Q - p2Q - p3Q));
  r.ray8 = (float)(-two34 * l->ray7 - two23 * l->ray11 + l->ray8 * (p3Q - p2Q - p4Q));
  r.ray9 =
      (float)(-two35 * l->ray7 - two24 * l->ray14 - two45 * l->ray8 - two23 * l->ray12 + l->ray9 * (p3Q + p4Q - p2Q));
  r.ray11 = (float)(two24 * l->ray7 - two23 * l->ray8 + l->ray11 * (p2Q - p3Q - p4Q));
  r.ray12 = (float)(two25 * l->ray7 - two23 * l->ray9 - two45 * l->ray11 + l->ray12 * (p2Q - p3Q + p4Q));
  r.ray14 =
      (float)(two35 * l->ray11 - two34 * l->ray12 - two24 * l->ray9 + two25 * l->ray8 + l->ray14 * (p2Q + p3Q - p4Q));

  r.normalize();
  return r;
}

float GA_CPU_Tracer::calculateDistanceBetweenPoints(Point* p0, Point* p1) {
  return sqrt((p0->point2 - p1->point2) * (p0->point2 - p1->point2) +
              (p0->point3 - p1->point3) * (p0->point3 - p1->point3) +
              (p0->point4 - p1->point4) * (p0->point4 - p1->point4));
}

Color GA_CPU_Tracer::traceRay(Ray* ray, Plane* rpl, Camera* cam, Color col, int recursion) {
  if (recursion == 0) {
    return col;
  }

  SplitScene* loc_scene = SplitScene::Instance();

  Point IP;
  Point actIP;
  double mindist = 10000.0f;
  bool intfound  = false;
  int objID      = 0;

  for (int j = 0; j < loc_scene->getMeshCnt(); j++) {
    // BoundingSphere Test
    if (!loc_scene->meshes[j].bound.doesRayIntersect2(ray))
      continue;

    for (int i = 0; i < loc_scene->meshes[j].getObjectCount(); i++) {
      Triangle act = loc_scene->meshes[j].objects[i];

      // Kein Schnitt
      if (!act.testIntersection(ray, rpl, actIP, true)) {
        continue;
      }

      intfound       = true;
      double actdist = calculateDistanceBetweenPoints(&actIP, ray->origin);

      if (actdist < mindist) {
        mindist = actdist;
        IP      = actIP;
        objID   = j;
      }
    }
  }

  if (!intfound)
    return col;

  Mesh2 actObj = loc_scene->meshes[objID];
  Material mat = loc_scene->mats[actObj.mat];

  // Licht berechnen
  Ray lightvec = Ray(&IP, &(loc_scene->light.pos));
  lightvec.normalize();
  Plane tpl           = Plane(&lightvec);
  Direction* lightdir = lightvec.getDirection();

  Color amb = calculateAmbient(&(loc_scene->light), &mat);
  col       = col + amb;

  Direction* eye = ray->getDirection();
  eye->invertDirection();

  // 3D
  Direction reflect = reflect3D(eye, &IP.n);
  reflect.normalize();
  Ray refray = Ray(&IP, &reflect);

  // 5D
  /*Ray refray = reflect5D(&lightvec, &IP.n);
  Direction* reflect = refray.getDirection();
  reflect->normalize();*/

  Plane refpl  = Plane(&refray);
  double NDotL = angleBetweenDirections(&IP.n, lightdir);
  if (NDotL > 0.0) {
    /*maxmesh=0;
    j=0;
    // Über Objekte iterieren
    for(int k = 0; k < loc_scene->objcnt; k++)
    {
      Sphere objsphere = loc_scene->objbounds[k];
      bool objindi = objsphere.doesRayIntersect2(&lightvec);

      // Objekt wir nicht geschnitten, also kann zum nächsten gesprungen werden
      this->testcnt++;
      if(!objindi)
      {
        j = j + loc_scene->objmeshes[k];
        maxmesh += loc_scene->objmeshes[k];
        continue;
      }

      maxmesh += loc_scene->objmeshes[k];
      Point tmp;
      for(; j<maxmesh; j++)
      {
        // BoundingSphere Test
        if(!loc_scene->meshes[j].bound.doesRayIntersect(&lightvec))
          continue;

        for (int i=0; i<loc_scene->meshes[j].getObjectCount(); i++)
        {
          Triangle act = loc_scene->meshes[j].objects[i];
          this->testcnt++;
          // Kein Schnitt
          if(act.testIntersection(&lightvec, &tpl, tmp, false))
          {
            col = traceRay(&refray, &refpl, cam, col, recursion-1);
            return col;
          }
        }
      }
    }*/

    // Shadowray auf Schnitt prüfen
    for (int j = 0; j < loc_scene->getMeshCnt(); j++) {
      // BoundingSphere Test
      if (!loc_scene->meshes[j].bound.doesRayIntersect(&lightvec)) {
        continue;
      }

      // Wenn in Sphere dann test gegen die Dreiecke davon
      Point tmp;
      for (int i = 0; i < loc_scene->meshes[j].getObjectCount(); i++) {
        // Wenn Schnitt
        if (loc_scene->meshes[j].objects[i].testIntersection(&lightvec, &tpl, tmp, false)) {
          col = traceRay(&refray, &refpl, cam, col, recursion - 1);
          return col;
        }
      }
    }

    Color dif = calculateDiffuse(&(loc_scene->light), &mat, (float)NDotL);
    col       = col + dif;

    double NdotHV = angleBetweenDirections(&reflect, lightdir);

    Color spec = calculateSpecular(&(loc_scene->light), &mat, (float)NdotHV);
    col        = col + spec;

    // delete reflect;
    delete eye;
  }

  // Rekursive Auführungen
  col = traceRay(&refray, &refpl, cam, col, recursion - 1);

  delete lightdir;
  return col;
}

void GA_CPU_Tracer::rayTraceSceneNew(void) {
  // Kamera instantiieren
  SplitScene* loc_scene = SplitScene::Instance();
  Camera* cam           = &(loc_scene->cam);

#pragma omp parallel for
  for (int x = 0; x < cam->imageWidth; x++) {
    for (int y = 0; y < cam->imageWidth; y++) {
      // ursprünglichen Sichtstrahl erstellen
      Ray ray   = cam->getCamRay(x, y);
      Plane rpl = Plane(&ray);

      Color col = traceRay(&ray, &rpl, cam, Color(0, 0, 0), 2);

      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 0] = (float)col.r;
      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 1] = (float)col.g;
      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 2] = (float)col.b;
    }
  }
}

int GA_CPU_Tracer::run() {
  return 0;
}

void GA_CPU_Tracer::renderImage() {
  if (rendered == false) {
    u_timeElapsed();
    rayTraceSceneNew();
    double time = u_timeElapsed();
    printf("Time needed for Rendering Scene: %f \n", time);
    Log::Instance()->Write("Time needed for Rendering Scene: %f", time);
    //		Log::Instance()->Write("Triangleintersections needed for Rendering this Scene: %d", this->testcnt);
    Log::Instance()->Write("Used Splitfac for this Scene: %d \n", SplitScene::Instance()->getSplitFac());
    rendered = true;
  }
}
