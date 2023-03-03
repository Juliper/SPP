#include "SplitScene.h"

#include "Transformation.h"
//#include "../Timing.h"
#include <omp.h>
#include <stdarg.h>

SplitScene* SplitScene::_instance = 0;

SplitScene* SplitScene::Instance() {
  if (_instance == 0) {
    _instance = new SplitScene();
  }
  return _instance;
}

SplitScene::SplitScene() {
  this->meshcnt      = 0;
  this->vtxcnt       = 0;
  this->openglvtxcnt = 0;
  this->matcount     = 0;
  this->oldmeshcnt   = 0;

  this->meshes.resize(0);
  this->vtx.resize(0);
  this->mats.resize(0);
  this->oldmeshes.resize(0);

  // Daten für Objekte
  this->objcnt = 0;
  this->objbounds.resize(0);
  this->objmeshes.resize(0);
  this->obj.resize(0);
  this->splitfac = 500;

  this->force_cpu_fallback = false;
}

int SplitScene::getMeshCnt() {
  return this->meshcnt;
}

int SplitScene::getOldMeshCnt() {
  return this->oldmeshcnt;
}

void SplitScene::addMesh(Mesh m) {
  this->oldmeshcnt++;
  m.configureMesh();
  this->oldmeshes.push_back(m);

  // Mesh muss Transformiert werden
  Mesh2 newmesh = Mesh2();

  // Parameter des alten kopieren
  newmesh.mat = m.mat;

  // Offset der Punkte auf die Faces addieren
  for (int i = 0; i < m.facecount; i++) {
    Face f = Face();
    f.p0   = m.faces[i].p0 + this->vtxcnt;
    f.p1   = m.faces[i].p1 + this->vtxcnt;
    f.p2   = m.faces[i].p2 + this->vtxcnt;

    newmesh.addFace(f);
  }

  // Punkte in das Array der Szene einlesen
  for (int i = 0; i < m.vtxcount; i++) {
    this->vtx.push_back(m.points[i]);
    this->vtxcnt++;
  }

  // Splitten und obj sichern
  //	u_timeElapsed();
  splitMesh(newmesh);
  //	double time = u_timeElapsed();
  //	printf("Time for loading scene from file: %f \n", time);
}

void SplitScene::configureScene() {
  openglvtx.resize(0);
  this->openglvtxcnt = this->vtxcnt;
  openglvtx.resize(this->openglvtxcnt);
  for (int i = 0; i < this->openglvtxcnt; i++) {
    openglvtx[i]   = Point(&(vtx[i]));
    openglvtx[i].n = Direction(vtx[i].n.direction2, vtx[i].n.direction3, vtx[i].n.direction4);
  }
  this->opengllight = this->light.pos;
}

void SplitScene::splitMesh(Mesh2 m) {
  // Daten für ganzes Objekt in der Szene speichern
  //#pragma omp parallel for
  for (int i = 0; i < m.facecount; i++) {
    m.addPoints(&(m.faces[i]));
  }

  // u_timeElapsed();

  // BoundingSphere für ganze Objekte
  BoundSphereCreator bsc = BoundSphereCreator();
  float* res             = bsc.createBound(&m, this->vtx);
  float r                = sqrt(res[3]);
  Point c                = Point(res[0], res[1], res[2]);

  Sphere b = Sphere(c, r);
  //	delete bsc;
  delete res;

  m.bound = b;
  this->obj.push_back(m);
  this->objbounds.push_back(b);
  this->objcnt++;

  int cnt            = 0;
  int meshes_for_obj = 0;
  Mesh2 tmpmesh;

  for (int i = 0; i < m.facecount; i++) {
    // Neues Netz Anfangen, da Maximalzahl der Faces im Schritt zuvor erreicht
    if (cnt == 0) {
      tmpmesh     = Mesh2();
      tmpmesh.mat = m.mat;
      cnt         = this->splitfac;
    }

    // Flächen kopieren
    tmpmesh.addFace(m.faces[i]);
    tmpmesh.addPoints(&(m.faces[i]));
    cnt--;

    if (cnt == 0) {
      tmpmesh.configureMesh(this->vtx);
      this->addConfiguredMesh(tmpmesh);
      // this->obj.push_back(tmpmesh);
      // this->objbounds.push_back(tmpmesh.bound);
      meshes_for_obj++;
    }
  }

  if (tmpmesh.facecount != 0) {
    tmpmesh.configureMesh(this->vtx);
    this->addConfiguredMesh(tmpmesh);
    // this->obj.push_back(tmpmesh);
    // this->objbounds.push_back(tmpmesh.bound);
    meshes_for_obj++;
  }

  this->objmeshes.push_back(meshes_for_obj);
  // this->objcnt++;
  // double time = u_timeElapsed();
  // printf("Time for configuring scene2: %f \n", time);
}

int SplitScene::getSplitFac() {
  return this->splitfac;
}

void SplitScene::setSplitFac(int fac) {
  this->splitfac = fac;
}

void SplitScene::addConfiguredMesh(Mesh2 m) {
  this->meshes.push_back(m);
  this->meshcnt++;
}

void SplitScene::addMaterial(Material mat) {
  this->mats.push_back(mat);
  this->matcount++;
}

int SplitScene::getMaterialCount() {
  return this->matcount;
}

void SplitScene::rotateScene(float angle, float x, float y, float z) {
  Transformation* rot = new Transformation();
  rot->setRotation(angle * 3.14f / 180.0f, x, y, z);

  for (int i = 0; i < this->vtxcnt; i++) {
    this->vtx[i] = rot->transformPoint(vtx[i]);
  }

  // Hier werden die BS der äußeren Kugeln passend umgelegt
  for (int i = 0; i < this->meshcnt; i++) {
    this->meshes[i].configureMesh(this->vtx);
  }

  // Für alle Objekte die Boundingsphere anpassen
  for (int i = 0; i < this->objcnt; i++) {
    this->obj[i].configureMesh(this->vtx);
    this->objbounds[i] = this->obj[i].bound;
  }

  // Licht ändern
  this->light.pos = rot->transformPoint(light.pos);

  // Boundingspheres ändern
  /*for(int i = 0; i<this->meshcnt; i++)
  {
    meshes[i].bound.MP = rot->transformPoint(meshes[i].bound.MP);
    meshes[i].bound.s2 = meshes[i].bound.MP.point2;
    meshes[i].bound.s3 = meshes[i].bound.MP.point3;
    meshes[i].bound.s4 = meshes[i].bound.MP.point4;
  }*/

  delete rot;
}

void SplitScene::resetScene() {
  for (int i = 0; i < this->vtxcnt; i++) {
    this->vtx[i] = openglvtx[i];
  }

  this->light.pos = this->opengllight;
}

void SplitScene::clearScene() {
  this->meshcnt      = 0;
  this->vtxcnt       = 0;
  this->openglvtxcnt = 0;
  this->matcount     = 0;
  this->oldmeshcnt   = 0;

  this->meshes.resize(0);
  this->vtx.resize(0);
  this->openglvtx.resize(0);
  this->mats.resize(0);
  this->oldmeshes.resize(0);

  // Daten für Objekte
  this->objcnt = 0;
  this->obj.resize(0);
  this->objbounds.resize(0);
  this->objmeshes.resize(0);
}

void SplitScene::TranslateSceneRTM() {
  rtm_stream.open("rtmscene.rtm", ios::out | ios::app);

  // Header für faces
  WriteRTMLine("number_of_faces %d", this->meshes[0].facecount);
  // Schreibe faces
  Mesh2 tmp = this->meshes[0];
  for (int i = 0; i < this->meshes[0].facecount; i++) {
    WriteRTMLine("face %d vtx %d %d %d", i, tmp.faces[i].p0, tmp.faces[i].p1, tmp.faces[i].p2);
  }

  WriteRTMLine("");

  // Header für vertices
  WriteRTMLine("number_of_vertices %d", this->meshes[0].vtxcount);
  // Schreibe vertices
  for (int i = 0; i < this->meshes[0].vtxcount; i++) {
    Point p = vtx[i];
    WriteRTMLine("vertex %i pos %4.10f %4.10f %4.10f normal %f %f %f uv 0 0", i, p.point2, p.point3, p.point4,
                 p.n.direction2, p.n.direction3, p.n.direction4);
  }

  // Farben und Co
  WriteRTMLine("");
  WriteRTMLine("ambient_color 0.19 0.07 0.02");
  WriteRTMLine("diffuse_color 0.7 0.27 0.08");
  WriteRTMLine("specular_color 0.25 0.13 0.08");
  WriteRTMLine("shininess 12.8");
  WriteRTMLine("");
  WriteRTMLine("direct_lighting 1.0");
  WriteRTMLine("reflection 1.0");
  WriteRTMLine("refraction 0.0 index 1.0");

  rtm_stream.close();
}

void SplitScene::exportTriangle(int i) {
  if (i > meshes[0].objectcount)
    return;

  Triangle act = meshes[0].objects[i];

  tria_stream.open("tria.txt", ios::out | ios::app);
  int temp = 0;
  temp     = (int)(act.p0.point2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [31:0] <= 32'h%08X;", temp);

  temp = (int)(act.p0.point3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [63:32] <= 32'h%08X;", temp);

  temp = (int)(act.p0.point4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [95:64] <= 32'h%08X;", temp);

  temp = (int)(act.p0.n.direction2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [127:96] <= 32'h%08X;", temp);

  temp = (int)(act.p0.n.direction3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [159:128] <= 32'h%08X;", temp);

  temp = (int)(act.p0.n.direction4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [191:160] <= 32'h%08X;", temp);

  // zweiter Punkt p1
  temp = (int)(act.p1.point2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [223:192] <= 32'h%08X;", temp);

  temp = (int)(act.p1.point3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [255:224] <= 32'h%08X;", temp);

  temp = (int)(act.p1.point4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [287:256] <= 32'h%08X;", temp);

  temp = (int)(act.p1.n.direction2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [319:288] <= 32'h%08X;", temp);

  temp = (int)(act.p1.n.direction3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [351:320] <= 32'h%08X;", temp);

  temp = (int)(act.p1.n.direction4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [383:352] <= 32'h%08X;", temp);

  // zweiter Punkt p2
  temp = (int)(act.p2.point2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [415:384] <= 32'h%08X;", temp);

  temp = (int)(act.p2.point3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [447:416] <= 32'h%08X;", temp);

  temp = (int)(act.p2.point4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [479:448] <= 32'h%08X;", temp);

  temp = (int)(act.p2.n.direction2 * SCALE);
  WriteTriaLine("TRIANGLE_IN [511:480] <= 32'h%08X;", temp);

  temp = (int)(act.p2.n.direction3 * SCALE);
  WriteTriaLine("TRIANGLE_IN [543:512] <= 32'h%08X;", temp);

  temp = (int)(act.p2.n.direction4 * SCALE);
  WriteTriaLine("TRIANGLE_IN [575:544] <= 32'h%08X;", temp);

  tria_stream.close();
}

void SplitScene::TranslateSceneMEM() {
  mem_stream.open("rtmscene.rtm", ios::out | ios::app);
  Mesh2 act_m = this->meshes[0];

  // Header für faces
  WriteMEMLine("00000004 %X", act_m.getObjectCount() * 18);

  for (int i = 0; i < act_m.getObjectCount(); i++) {
    Triangle act_t = act_m.objects[i];
    int temp       = 0;
    // Dreieck nach und nach in Datei schreiben:
    // Erster Punkt p0
    temp           = (int)(act_t.p0.point2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p0.point3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p0.point4 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p0.n.direction2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p0.n.direction3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p0.n.direction4 * SCALE);
    Log::Instance()->Write("%08X", temp);

    // zweiter Punkt p1
    temp = (int)(act_t.p1.point2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p1.point3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p1.point4 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p1.n.direction2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p1.n.direction3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p1.n.direction4 * SCALE);
    Log::Instance()->Write("%08X", temp);

    // zweiter Punkt p2
    temp = (int)(act_t.p2.point2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p2.point3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p2.point4 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p2.n.direction2 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p2.n.direction3 * SCALE);
    Log::Instance()->Write("%08X", temp);

    temp = (int)(act_t.p2.n.direction4 * SCALE);
    Log::Instance()->Write("%08X", temp);
  }

  rtm_stream.close();
}

void SplitScene::WriteRTMLine(const char* logline, ...) {
  va_list argList;
  char cbuffer[1024];
  va_start(argList, logline);
  vsnprintf(cbuffer, 1024, logline, argList);
  va_end(argList);
  rtm_stream << cbuffer << endl;
}

void SplitScene::WriteMEMLine(const char* logline, ...) {
  va_list argList;
  char cbuffer[1024];
  va_start(argList, logline);
  vsnprintf(cbuffer, 1024, logline, argList);
  va_end(argList);
  mem_stream << cbuffer << endl;
}

void SplitScene::WriteTriaLine(const char* logline, ...) {
  va_list argList;
  char cbuffer[1024];
  va_start(argList, logline);
  vsnprintf(cbuffer, 1024, logline, argList);
  va_end(argList);
  tria_stream << cbuffer << endl;
}