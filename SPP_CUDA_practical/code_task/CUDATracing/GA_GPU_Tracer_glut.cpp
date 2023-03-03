#include "GA_GPU_Tracer_glut.h"

#include "MeshSphere.h"
#include "Transformation.h"

#include <math.h>
//#include <GL/glut.h>
#include "Log.h"
#include "kernel.cuh"

#include <fstream>
#include <iostream>
#include <string>

std::string GA_GPU_Tracer::convertToString(const char* filename) {
  size_t size;
  char* str;
  std::string s;

  std::fstream f(filename, (std::fstream::in | std::fstream::binary));

  if (f.is_open()) {
    size_t fileSize;
    f.seekg(0, std::fstream::end);
    size = fileSize = (size_t)f.tellg();
    f.seekg(0, std::fstream::beg);

    str = new char[size + 1];
    if (!str) {
      f.close();
      return NULL;
    }

    f.read(str, fileSize);
    f.close();
    str[size] = '\0';

    s = str;

    return s;
  }
  return "";
}

int GA_GPU_Tracer::cleanupCL(void) {
  /*cl_int status;

  status = clReleaseKernel(kernel);
  if(status != CL_SUCCESS)
{
  std::cout<<"Error: In clReleaseKernel \n";
  return 1;
}
  status = clReleaseProgram(program);
  if(status != CL_SUCCESS)
{
  std::cout<<"Error: In clReleaseProgram\n";
  return 1;
}

  status = clReleaseCommandQueue(commandQueue);
  if(status != CL_SUCCESS)
{
  std::cout<<"Error: In clReleaseCommandQueue\n";
  return 1;
}
  status = clReleaseContext(context);
  if(status != CL_SUCCESS)
{
  std::cout<<"Error: In clReleaseContext\n";
  return 1;
}*/

  return 0;
}

// virtuelle Methoden

// Used to build the scene
int GA_GPU_Tracer::setup() {
  this->scene = SplitScene::Instance();
  // Variablenspeicherbedarf ermitteln
  facecnt     = 0;
  vtxcnt      = 0;

  for (int i = 0; i < scene->getMeshCnt(); i++) {
    // actual mesh
    facecnt += scene->meshes[i].facecount;
  }

  // Speicherplätze für Variablen, die beschrieben werden instantiierten
  faces         = (int*)malloc(facecnt * 4 * sizeof(int));
  vertices      = (float*)malloc(this->scene->vtxcnt * 4 * sizeof(float));
  vertices_norm = (float*)malloc(this->scene->vtxcnt * 4 * sizeof(float));
  loc_counts    = (int*)malloc(this->scene->getMeshCnt() * sizeof(int));
  spheres       = (float*)malloc(this->scene->getMeshCnt() * 4 * sizeof(float));
  materials     = (float*)malloc(this->scene->getMaterialCount() * 16 * sizeof(float));
  matidx        = (int*)malloc(this->scene->getMeshCnt() * sizeof(int));

  objmeshes = (int*)malloc(this->scene->objcnt * sizeof(int));
  objbounds = (float*)malloc(this->scene->objcnt * 4 * sizeof(float));

  // Copy Values to Bufferarrays
  int absfcnt = 0;

  for (int i = 0; i < scene->objcnt; i++) {
    objmeshes[i] = scene->objmeshes[i];

    // LA
    if (this->type == 1) {
      objbounds[4 * i + 0] = scene->objbounds[i].MP.point2;
      objbounds[4 * i + 1] = scene->objbounds[i].MP.point3;
      objbounds[4 * i + 2] = scene->objbounds[i].MP.point4;
      objbounds[4 * i + 3] = scene->objbounds[i].radius;
    }

    // GA
    else {
      objbounds[4 * i + 0] = scene->objbounds[i].s2;
      objbounds[4 * i + 1] = scene->objbounds[i].s3;
      objbounds[4 * i + 2] = scene->objbounds[i].s4;
      objbounds[4 * i + 3] = scene->objbounds[i].s5;
    }
  }

  for (int i = 0; i < scene->getMaterialCount(); i++) {
    // Materials
    // ambient
    materials[16 * i + 3]  = scene->mats[i].ambient.r;
    materials[16 * i + 4]  = scene->mats[i].ambient.g;
    materials[16 * i + 5]  = scene->mats[i].ambient.b;
    // diffuse
    materials[16 * i + 6]  = scene->mats[i].diffuse.r;
    materials[16 * i + 7]  = scene->mats[i].diffuse.g;
    materials[16 * i + 8]  = scene->mats[i].diffuse.b;
    // specular
    materials[16 * i + 9]  = scene->mats[i].specular.r;
    materials[16 * i + 10] = scene->mats[i].specular.g;
    materials[16 * i + 11] = scene->mats[i].specular.b;
    // shininess
    materials[16 * i + 12] = scene->mats[i].shininess;
  }

  for (int i = 0; i < scene->getMeshCnt(); i++) {
    // Read actual Mesh
    Mesh2 act_mesh = scene->meshes[i];

    // Counts
    loc_counts[i] = act_mesh.facecount;

    // Bounding Spheres
    // LA
    if (this->type == 1) {
      spheres[4 * i + 0] = act_mesh.bound.MP.point2;
      spheres[4 * i + 1] = act_mesh.bound.MP.point3;
      spheres[4 * i + 2] = act_mesh.bound.MP.point4;
      spheres[4 * i + 3] = act_mesh.bound.radius;
    }

    // GA
    else {
      spheres[4 * i + 0] = act_mesh.bound.s2;
      spheres[4 * i + 1] = act_mesh.bound.s3;
      spheres[4 * i + 2] = act_mesh.bound.s4;
      spheres[4 * i + 3] = act_mesh.bound.s5;
    }

    // Material-Indices
    matidx[i] = act_mesh.mat;

    // Faces
    for (int j = 0; j < act_mesh.facecount; j++) {
      faces[4 * absfcnt + 0] = act_mesh.faces[j].p0;
      faces[4 * absfcnt + 1] = act_mesh.faces[j].p1;
      faces[4 * absfcnt + 2] = act_mesh.faces[j].p2;
      faces[4 * absfcnt + 3] = 0;  // dummy
      absfcnt++;
    }
  }

  for (int i = 0; i < this->scene->vtxcnt; i++) {
    // P from List
    vertices[4 * i + 0] = this->scene->vtx[i].point2;
    vertices[4 * i + 1] = this->scene->vtx[i].point3;
    vertices[4 * i + 2] = this->scene->vtx[i].point4;
    vertices[4 * i + 3] = 0.0f;  // dummy

    // Normal of this point
    vertices_norm[4 * i + 0] = this->scene->vtx[i].n.direction2;
    vertices_norm[4 * i + 1] = this->scene->vtx[i].n.direction3;
    vertices_norm[4 * i + 2] = this->scene->vtx[i].n.direction4;
    vertices_norm[4 * i + 3] = 0.0f;  // dummy
  }

  return 0;
}

int GA_GPU_Tracer::verifyResults() {
  return 1;
}

int GA_GPU_Tracer::initialize() {
  return 0;
}

void GA_GPU_Tracer::printStats() {
  return;
}

int GA_GPU_Tracer::cleanup() {
  return 1;
}

int GA_GPU_Tracer::setupCL() {
  return 0;
}

int GA_GPU_Tracer::runCLKernels() {
  return 0;
}

int GA_GPU_Tracer::run() {
  return 0;
}

void GA_GPU_Tracer::renderImage() {
  if (!rendered) {
    setup();

    float* result = (float*)malloc(width * height * 4 * sizeof(float));

    traceWithCuda(result, width, height, scene, width, height, facecnt, spheres, faces, vertices, vertices_norm,
                  loc_counts, materials, matidx, objbounds, objmeshes);

    // Copy Colorvalues to Bitmap
    int cnt = 0;
    for (int x = 0; x < width; x++) {
      for (int y = 0; y < height; y++) {
        float r = result[(y * (int)width + x) * 4 + 0];
        if (r < 0.0f)
          r = 0.0f;
        if (r > 1.0f)
          r = 1.0f;

        float g = result[(y * (int)width + x) * 4 + 1];
        if (g < 0.0f)
          g = 0.0f;
        if (g > 1.0f)
          g = 1.0f;

        float b = result[(y * (int)width + x) * 4 + 2];
        if (b < 0.0f)
          b = 0.0f;
        if (b > 1.0f)
          b = 1.0f;

        checkImage[((height - y - 1) * width + x) * 3 + 0] = r;
        checkImage[((height - y - 1) * width + x) * 3 + 1] = g;
        checkImage[((height - y - 1) * width + x) * 3 + 2] = b;

        cnt += (int)result[(y * (int)width + x) * 4 + 3];
      }
    }
    free(result);
    /*setupCL();
    runCLKernels();
    cleanupCL();*/
    rendered = true;
  }
  return;
}
