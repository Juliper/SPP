
//#define USE_GLUT
#ifdef USE_GLUT
#include <GL/glut.h>
#endif

#include "BitmapWriter.h"
#include "Camera.h"
#include "Color.h"
#include "Constants.h"
#include "DirectionPoint.h"
#include "GA_CPU_Tracer_glut.h"
#include "Material.h"
#include "MeshFileReader.h"
#include "SplitScene.h"
#include "Transformation.h"
//#include "GA_GPU_Tracer_glut.h"
#include "BoundSphereCreator.h"
#include "LA_CPU_Tracer_glut.h"
#include "LA_GPU_Tracer_glut.h"
#include "Log.h"
#include "MaterialBuilder.h"
#include "timing.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

GA_CPU_Tracer* gacputr;
// GA_GPU_Tracer* gagputr;
LA_CPU_Tracer* lacputr;
LA_GPU_Tracer* lagputr;

int tracer_type = 0;

void instantiateScene(SplitScene* scene, const char* filename, int splitfac) {
  MaterialBuilder* mb    = MaterialBuilder::Instance();
  MeshFileReader* reader = new MeshFileReader();
  Mesh mesh0;

  scene->setSplitFac(splitfac);

  fstream fin(filename);
  string buffer;
  string tmp_array[200];
  int i      = 0;
  float rotZ = 0.0f;
  float rotY = 0.0f;
  float rotX = 0.0f;

  while (fin.good()) {
    getline(fin, buffer);

    // Wenn es ein Mesh ist
    if (buffer == "MESH") {
      u_timeElapsed();
      getline(fin, buffer);
      float transx = (float)atof(buffer.c_str());
      getline(fin, buffer);
      float transy = (float)atof(buffer.c_str());
      getline(fin, buffer);
      float transz = (float)atof(buffer.c_str());

      Transformation* trans0 = new Transformation();
      trans0->setTranslation(transx, transy, transz);

      getline(fin, buffer);
      char* filep = new char[1024];
      std::strcpy(filep, buffer.c_str());

      getline(fin, buffer);
      int ttpye = atoi(buffer.c_str());

      // Material lesen
      getline(fin, buffer);
      int matidx = atoi(buffer.c_str());

      // Szene konfigurieren
      mesh0 = reader->readMeshFromFile(filep, trans0, ttpye);

      int matcnt = scene->getMaterialCount();
      mesh0.mat  = matcnt;
      scene->addMaterial(MaterialBuilder::Instance()->getMaterial(matidx));

      scene->addMesh(mesh0);

    }

    else if (buffer == "ROTATION") {
      getline(fin, buffer);
      rotZ = (float)atof(buffer.c_str());
      getline(fin, buffer);
      rotY = (float)atof(buffer.c_str());
      getline(fin, buffer);
      rotX = (float)atof(buffer.c_str());
    }
  }

  fin.close();

  // Configure Camera
  Point* campos = new Point(0.0f, 0.0f, 4.0f);
  scene->cam    = Camera(campos, 0.74f, 0.74f, IMAGEWIDTH, IMAGEHEIGHT);

  // Configure Light
  Point lightpos = Point(0.0f, 0.0f, 4.0f);
  scene->light = Light(lightpos, new Color(1.0f, 1.0f, 1.0f), new Color(1.0f, 1.0f, 1.0f), new Color(1.0f, 1.0f, 1.0f));
  scene->configureScene();

  // Szene rotieren
  scene->resetScene();
  scene->rotateScene(rotZ, 0, 0, 1);
  scene->rotateScene(rotY, 0, 1, 0);
  scene->rotateScene(rotX, 1, 0, 0);
  return;
}

#ifdef USE_GLUT
void init(void) {
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glShadeModel(GL_SMOOTH);
}

void reshape(int w, int h) {
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
}

void display(void) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glRasterPos2i(0, 0);

  if (tracer_type == 1 || tracer_type == 5) {
    // glDrawPixels(IMAGEWIDTH, IMAGEHEIGHT, GL_RGB,
    // GL_FLOAT, gagputr->checkImage);
  } else if (tracer_type == 2 || tracer_type == 6)
    glDrawPixels(IMAGEWIDTH, IMAGEHEIGHT, GL_RGB, GL_FLOAT, lagputr->checkImage);
  else if (tracer_type == 3)
    glDrawPixels(IMAGEWIDTH, IMAGEHEIGHT, GL_RGB, GL_FLOAT, gacputr->checkImage);
  else if (tracer_type == 4)
    glDrawPixels(IMAGEWIDTH, IMAGEHEIGHT, GL_RGB, GL_FLOAT, lacputr->checkImage);
  else
    glDrawPixels(IMAGEWIDTH, IMAGEHEIGHT, GL_RGB, GL_FLOAT, lacputr->checkImage);

  glFlush();
}
#endif

int main(int argc, char* argv[]) {
  int obj      = 1;
  tracer_type  = 2;
  int show_pic = 0;
  const char* filename;
  std::string sfilename;

  printf("Es wurden %d Parameter angegeben.\n\n", argc);

  if (argc < 4) {
    printf("Falsche Anzahl von Parametern! Tracerdurchlauf abgebrochen.");
    exit(0);
  }

  // Parameter einlesen
  for (int i = 0; i < argc; i++) {
    if (i == 1)
      tracer_type = atoi(argv[1]);
    if (i == 2) {
      filename  = argv[2];
      sfilename = std::string(filename);
      obj       = atoi(argv[2]);
    }
    if (i == 3)
      show_pic = atoi(argv[3]);
  }

  // Splitfaktor auf default lassen oder einlesen, falls vom Benutzer übergeben
  int splitfac = 500;
  if (argc >= 5) {
    splitfac = atoi(argv[4]);
  }

  printf("Splitfac ist %d.\n", splitfac);

  string begin("LOG--");
  begin += sfilename;
  Log::Instance(begin)->Write("Starting Visual C++-Raytracer of type %i with SceneFile: %s.", tracer_type, filename);
  instantiateScene(SplitScene::Instance(), filename, splitfac);

  if (tracer_type == 1 || tracer_type == 2) {
    SplitScene::Instance()->force_cpu_fallback = false;
    lagputr                                    = new LA_GPU_Tracer(IMAGEWIDTH, IMAGEHEIGHT);
    lagputr->renderImage();
    writeBMPToFile(lagputr->checkImage, IMAGEWIDTH, IMAGEHEIGHT);
  }

  else if (tracer_type == 3) {
    gacputr = new GA_CPU_Tracer(IMAGEWIDTH, IMAGEHEIGHT);
    gacputr->renderImage();
    writeBMPToFile(gacputr->checkImage, IMAGEWIDTH, IMAGEHEIGHT);
  }

  else if (tracer_type == 4) {
    lacputr = new LA_CPU_Tracer(IMAGEWIDTH, IMAGEHEIGHT);
    lacputr->renderImage();
    writeBMPToFile(lacputr->checkImage, IMAGEWIDTH, IMAGEHEIGHT);
  }

  else {
    lacputr = new LA_CPU_Tracer(IMAGEWIDTH, IMAGEHEIGHT);
    lacputr->renderImage();
    writeBMPToFile(lacputr->checkImage, IMAGEWIDTH, IMAGEHEIGHT);
  }

#ifdef USE_GLUT
  if (show_pic) {
    // Draw calculated Image
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(IMAGEWIDTH, IMAGEHEIGHT);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("CPU Tracer with use of Intel Compiler");
    init();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();

    return 0;
  }
#endif
}
