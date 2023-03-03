#include "LA_CPU_Tracer_glut.h"

#include "Camera.h"
#include "Color.h"
#include "Colorfunctions.h"
#include "DirectionPoint.h"
#include "Light.h"
#include "Log.h"
#include "MeshSphere.h"
#include "float.h"
#include "timing.h"

#include <math.h>
#include <omp.h>
//#include <GL/glut.h>

LA_CPU_Tracer::LA_CPU_Tracer(int width, int height) : virtualTracer(width, height) {
  this->width      = width;
  this->height     = height;
  this->checkImage = new float[width * height * 3];
}

Direction LA_CPU_Tracer::reflect3DLA(Direction* l, Direction* n) {
  l->normalize();
  n->normalize();

  float s = (*l) % (*n);

  return Direction(2.0f * n->direction2 * s - l->direction2, 2.0f * n->direction3 * s - l->direction3,
                   2.0f * n->direction4 * s - l->direction4);
}

float LA_CPU_Tracer::calculateDistanceBetweenPoints3D(Point* p0, Point* p1) {
  return sqrt((p1->point2 - p0->point2) * (p1->point2 - p0->point2) +
              (p1->point3 - p0->point3) * (p1->point3 - p0->point3) +
              (p1->point4 - p0->point4) * (p1->point4 - p0->point4));
}

Color LA_CPU_Tracer::traceRay(Ray3D* ray, Camera* cam, Color col, int recu) {
  if (recu == 0) {
    return col;
  }

  SplitScene* loc_scene = SplitScene::Instance();

  Point IP;
  Point actIP;
  float mindist = FLT_MAX;
  bool intfound = false;
  int objID     = 0;

  for (int j = 0; j < loc_scene->getMeshCnt(); j++) {
    // BoundingSphere Test
    // testcnt++;
    if (!loc_scene->meshes[j].bound.doesRayIntersect(ray))
      continue;

    // Wenn in Sphere dann test gegen die Dreiecke davon
    for (int i = 0; i < loc_scene->meshes[j].getObjectCount(); i++) {
      Triangle act = loc_scene->meshes[j].objects[i];

      // Kein Schnitt
      if (!act.testIntersection(ray)) {
        continue;
      }

      // Wenn Schnitt prüfen ob der näheste
      intfound      = true;
      actIP         = act.findIntersectionPoint();
      float actdist = calculateDistanceBetweenPoints3D(&actIP, &(ray->ori));

      if (actdist < mindist) {
        mindist = actdist;
        IP      = actIP;
        objID   = j;
      }
    }
  }

  if (!intfound) {
    return col;
  }

  /*Point IP;
  Point actIP;
  float mindist = FLT_MAX;
  bool intfound = false;
  int objID = 0;
  int j=0;
  int maxmesh=0;

  // Über Objekte iterieren
  for(int k = 0; k < loc_scene->objcnt; k++)
  {
    Sphere objsphere = loc_scene->objbounds[k];
    bool objindi = objsphere.doesRayIntersect(ray);

    // Objekt wir nicht geschnitten, also kann zum nächsten gesprungen werden
    if(!objindi)
    {
      j = j + loc_scene->objmeshes[k];
      maxmesh += loc_scene->objmeshes[k];
      continue;
    }

    maxmesh += loc_scene->objmeshes[k];
    for(; j<maxmesh; j++)
    {
      // BoundingSphere Test
      if(!loc_scene->meshes[j].bound.doesRayIntersect(ray))
        continue;

      for (int i=0; i<loc_scene->meshes[j].getObjectCount(); i++)
      {
        Triangle act = loc_scene->meshes[j].objects[i];

        // Kein Schnitt
        if(!act.testIntersection(ray))
        {
          continue;
        }

        intfound=true;
        actIP = act.findIntersectionPoint();
        float actdist = calculateDistanceBetweenPoints3D(&actIP, &(ray->ori));

        if(actdist < mindist)
        {
          mindist = actdist;
          IP = actIP;
          objID = j;
        }
      }
    }
  }

  if(!intfound)
    return col;*/

  Direction eye = ray->dir;
  eye.invertDirection();
  Direction eyeref = reflect3DLA(&eye, &IP.n);
  eyeref.normalize();
  Ray3D refray = Ray3D(IP, eyeref);

  // Licht berechnen
  Ray3D lightvec     = Ray3D(IP, loc_scene->light.pos);
  Direction lightdir = lightvec.dir;

  Mesh2 actObj = loc_scene->meshes[objID];
  Material mat = loc_scene->mats[actObj.mat];

  Color amb = calculateAmbient(&(loc_scene->light), &mat);
  col       = col + amb;

  float NDotL = angleBetweenDirections(&IP.n, &lightdir);
  if (NDotL > 0.0) {
    /*maxmesh=0;
    j=0;
    // Über Objekte iterieren
    for(int k = 0; k < loc_scene->objcnt; k++)
    {
      Sphere objsphere = loc_scene->objbounds[k];
      bool objindi = objsphere.doesRayIntersect(&lightvec);

      // Objekt wir nicht geschnitten, also kann zum nächsten gesprungen werden
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
          // Kein Schnitt
          if(loc_scene->meshes[j].objects[i].testIntersection(&lightvec))
          {
            col = traceRay(&refray,cam,col,recu-1);;
            return col;
          }
        }
      }
    }*/

    // Shadowray auf Schnitt prüfen
    for (int j = 0; j < loc_scene->getMeshCnt(); j++) {
      // BoundingSphere Test
      // testcnt++;
      if (!loc_scene->meshes[j].bound.doesRayIntersect(&lightvec))
        continue;

      // Wenn in Sphere dann test gegen die Dreiecke davon
      for (int i = 0; i < loc_scene->meshes[j].getObjectCount(); i++) {
        // Wenn Schnitt
        if (loc_scene->meshes[j].objects[i].testIntersection(&lightvec)) {
          col = traceRay(&refray, cam, col, recu - 1);
          return col;
        }
      }
    }

    Color dif = calculateDiffuse(&(loc_scene->light), &mat, NDotL);
    col       = col + dif;

    float NdotHV = angleBetweenDirections(&eyeref, &lightdir);

    Color spec = calculateSpecular(&(loc_scene->light), &mat, NdotHV);
    col        = col + spec;
  }

  // recursive calls
  col = traceRay(&refray, cam, col, recu - 1);
  return col;
}

void LA_CPU_Tracer::rayTraceScene(void) {
  SplitScene* loc_scene = SplitScene::Instance();
  // Kamera instantiieren
  Camera* cam           = &(loc_scene->cam);

#pragma omp parallel for
  for (int x = 0; x < cam->imageWidth; x++) {
    for (int y = 0; y < cam->imageHeight; y++) {
      // ursprünglichen Sichtstrahl erstellen
      Ray3D ray = cam->getCamRay3D(x, y);
      Color col = traceRay(&ray, cam, Color(0, 0, 0), 2);

      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 0] = (float)col.r;
      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 1] = (float)col.g;
      checkImage[((cam->imageWidth - y - 1) * cam->imageHeight + x) * 3 + 2] = (float)col.b;
    }
  }
}

int LA_CPU_Tracer::run() {
  return 0;
}

void LA_CPU_Tracer::renderImage() {
  if (rendered == false) {
    u_timeElapsed();
    rayTraceScene();
    double time = u_timeElapsed();
    printf("Time needed for Rendering Scene: %f \n", time);
    Log::Instance()->Write("Time needed for Rendering Scene: %f \n", time);
    rendered = true;
  }
}
