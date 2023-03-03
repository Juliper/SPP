#include "kernel.cuh"

#include <iostream>
#include <stdio.h>

// Constant buffers
// TODO Aufgabe 3)
__constant__ float camera_Buffer[8];
__constant__ float lightdat_Buffer[12];

// Normal buffers
float4* result_Buffer;
int2 resinfo_Buffer;
float4* spheres_Buffer;
int4* faces_Buffer;
float4* vertices_Buffer;
float4* normals_Buffer;
int* counts_Buffer;
float* materials_Buffer;
int* matidx_Buffer;
float4* objbounds_Buffer;
int* objmeshes_Buffer;
int4* result_AnaVals;

cudaError_t setupCuda() {
  cudaError_t cudaStatus;
  cudaStatus = cudaSetDevice(0);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
  }
  return cudaStatus;
}

////////// END BLOCK: DO NOT TOUCH ANY CODE IN THIS HELPER FUNCTIONS

/*
 * width = Horizontal no. of pixels of rendered image
 * height = Vertical no. of pixels of rendered image
 */
cudaError_t setupCudaMemory(SplitScene* scene, const int width, const int height, const int facecnt, float* spheres,
                            int* faces, float* vertices, float* normals, int* loc_counts, float* materials, int* matidx,
                            float* objbounds, int* objmeshes) {
  cudaError_t cudaStatus;
  size_t size;

  // Allocate and fill GPU buffers (faces, vertices, result, ...)    .
  // Buffer for result picture

  // TODO Aufgabe 1)
  cudaStatus = cudaMalloc((void**)&result_Buffer, width * height * sizeof(float4));
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "setupCudaMemory failed! Couldn't initialize \"result_Buffer\".");
    return cudaStatus;
  }
  cudaMemset((void*)result_Buffer, 0, width * height * sizeof(float4));
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "setupMemset failed! Couldn't initialize \"result_Buffer\" with 0.");
    return cudaStatus;
  }

  // Buffer for vertices
  // TODO Aufgabe 2)
  cudaStatus = cudaMalloc((void**)&vertices_Buffer, scene->vtxcnt * sizeof(float4));
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "setupCudaMemory failed! Couldn't initialize \"vertices_Buffer\".");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(vertices_Buffer, vertices, scene->vtxcnt * sizeof(float4), cudaMemcpyHostToDevice);
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "cudaMemcpy failed! Couldn't copy \"vertices\" to \"vertices_Buffer\" from Host to Device");
    return cudaStatus;
  }

  // Buffer for camera
  float camera[8];
  camera[0] = scene->cam.pos->point2;      // x
  camera[1] = scene->cam.pos->point3;      // y
  camera[2] = scene->cam.pos->point4;      // z
  camera[3] = scene->cam.cameraFOVwidth;   // cameraFOVwidth
  camera[4] = scene->cam.cameraFOVheight;  // cameraFOVheight
  camera[5] = scene->cam.pixelWidth;       // pixelWidth
  camera[6] = scene->cam.pixelHeight;      // camera.s[4] / height; // pixelHeight
  camera[7] = 0.0f;                        // dummy

  // Copy camera to constant camera_Buffer
  // TODO Aufgabe 3)
  cudaStatus = cudaMemcpyToSymbol(camera_Buffer, &camera, 8 * sizeof(float));
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "cudaMemcpy failed! Couldn't copy \"camera\" to \"camera_Buffer\" from Host to Device");
    return cudaStatus;
  }

  // Buffer for light
  float lightdat[12];
  lightdat[0]  = scene->light.pos.point2;
  lightdat[1]  = scene->light.pos.point3;
  lightdat[2]  = scene->light.pos.point4;
  lightdat[3]  = scene->light.ambient->r;
  lightdat[4]  = scene->light.ambient->g;
  lightdat[5]  = scene->light.ambient->b;
  lightdat[6]  = scene->light.diffuse->r;
  lightdat[7]  = scene->light.diffuse->g;
  lightdat[8]  = scene->light.diffuse->b;
  lightdat[9]  = scene->light.specular->r;
  lightdat[10] = scene->light.specular->g;
  lightdat[11] = scene->light.specular->b;

  // Copy light to constant light_Buffer

  // TODO Aufgabe 3)
  cudaStatus = cudaMemcpyToSymbol(lightdat_Buffer, &lightdat, 12 * sizeof(float));
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "cudaMemcpy failed! Couldn't copy \"lightdat\" to \"lightdat_Buffer\" from Host to Device");
    return cudaStatus;
  }

  /////// DO NOT TOUCH CODE OF THIS FUNCTION FROM HERE ON! ///////
  // Buffer for resolution info
  resinfo_Buffer.x = width;
  resinfo_Buffer.y = height;

  // Buffer for bounding spheres
  size       = scene->getMeshCnt();
  cudaStatus = cudaMalloc((void**)&spheres_Buffer, size * sizeof(float4));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for spheres_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(spheres_Buffer, spheres, size * sizeof(float) * 4, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for spheres_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for faces
  size       = facecnt;
  cudaStatus = cudaMalloc((void**)&faces_Buffer, size * sizeof(int4));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for faces_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(faces_Buffer, faces, size * sizeof(int) * 4, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for faces_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for normals
  size       = scene->vtxcnt;
  cudaStatus = cudaMalloc((void**)&normals_Buffer, size * sizeof(float4));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for normals_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(normals_Buffer, normals, size * sizeof(float) * 4, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for normals_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for counts
  size       = scene->getMeshCnt();
  cudaStatus = cudaMalloc((void**)&counts_Buffer, size * sizeof(int));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for counts_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(counts_Buffer, loc_counts, size * sizeof(int), cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for counts_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for materials
  size       = scene->getMaterialCount() * 16;
  cudaStatus = cudaMalloc((void**)&materials_Buffer, size * sizeof(float));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for materials_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(materials_Buffer, materials, size * sizeof(float), cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for materials_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for material indices
  size       = scene->getMeshCnt();
  cudaStatus = cudaMalloc((void**)&matidx_Buffer, size * sizeof(int));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for matidx_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(matidx_Buffer, matidx, size * sizeof(int), cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for matidx_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for object bounds
  size       = scene->objcnt;
  cudaStatus = cudaMalloc((void**)&objbounds_Buffer, size * sizeof(float4));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for objbounds_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(objbounds_Buffer, objbounds, size * sizeof(float) * 4, cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for objbounds_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for object meshes
  size       = scene->objcnt;
  cudaStatus = cudaMalloc((void**)&objmeshes_Buffer, size * sizeof(int));
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMalloc for objmeshes_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaMemcpy(objmeshes_Buffer, objmeshes, size * sizeof(int), cudaMemcpyHostToDevice);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy for objmeshes_Buffer failed!");
    return cudaStatus;
  }

  return cudaStatus;
}

cudaError_t freeCudaMemory() {
  cudaError_t cudaStatus;

  // Free GPU buffers (faces, vertices, result, ...)
  // Buffer for result picture

  // TODO Aufgabe 9
  cudaStatus = cudaFree(result_Buffer);
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "freeCUDA for result_Buffer failed!");
    return cudaStatus;
  }
  cudaStatus = cudaFree(vertices_Buffer);
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "freeCUDA for vertices_Buffer failed!");
    return cudaStatus;
  }

  /////// DO NOT TOUCH FUNCTION CODE FROM HERE ON! ///////

  // Buffer for bounding spheres
  cudaStatus = cudaFree(spheres_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for spheres_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for faces
  cudaStatus = cudaFree(faces_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for faces_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for normals
  cudaStatus = cudaFree(normals_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for normals_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for counts
  cudaStatus = cudaFree(counts_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for counts_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for materials
  cudaStatus = cudaFree(materials_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for materials_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for material indices
  cudaStatus = cudaFree(matidx_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for matidx_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for object bounds
  cudaStatus = cudaFree(objbounds_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for objbounds_Buffer failed!");
    return cudaStatus;
  }

  // Buffer for object meshes
  cudaStatus = cudaFree(objmeshes_Buffer);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCUDA for objmeshes_Buffer failed!");
    return cudaStatus;
  }

  return cudaStatus;
}

__device__ float4 normalizeDirection(float4 d) {
  float4 nd;
  // TODO Aufgabe 6
  float vectorAmount = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
  nd.x = (1 / vectorAmount) * d.x;
  nd.y = (1 / vectorAmount) * d.y;
  nd.z = (1 / vectorAmount) * d.z;

  d.w = 0.0f;

  return nd;
}

__device__ float4 invertDirection(float4 d) {
  float4 nd;
  nd.x = -d.x;
  nd.y = -d.y;
  nd.z = -d.z;
  return nd;
}

__device__ float distanceBetweenPoints(float4 p0, float4 p1) {
  return sqrtf((p0.x - p1.x) * (p0.x - p1.x) + (p0.y - p1.y) * (p0.y - p1.y) + (p0.z - p1.z) * (p0.z - p1.z));
}

__device__ float4 getCamRay3D(const float* cam, int x, int y, int2 resinfo) {
  float4 d;
  d.x = x * cam[5] - 0.5f * resinfo.x * cam[5];
  d.y = y * cam[6] - 0.5f * resinfo.y * cam[6];
  d.z = -1.0f;
  d.w = 0.0f;
  d   = normalizeDirection(d);
  return d;
}

__device__ float4 reflect3DLA(float4 l, float4 n) {
  float4 r;
  // TODO Aufgabe 7
  float scalar = l.x * n.x + l.y * n.y + l.z * n.z;
  r.x = 2 * scalar * n.x - l.x;
  r.y = 2 * scalar * n.y - l.y;
  r.z = 2 * scalar * n.z - l.z;

  r.w = 0.0f;

  return r;
}

/*
Tests whether a 3D Ray hits a bounding Sphere
*/
__device__ bool rayIntersectsSphere(float4 s, float4 r, float4 o) {
  // TODO Aufgabe 8

  // Ist die Diskriminante ß^2 - y nichtnegativ (-> mind. eine Lösung) oder
  // negativ (-> keine Lösung)
  float radius = s.w; 
  float4 diff;    // Wert von o - c
  diff.x = o.x - s.x;
  diff.y = o.y - s.y;
  diff.z = o.z - s.z;
  float beta = r.x * diff.x + r.y * diff.y + r.z * diff.z;
  float gamma = (diff.x * diff.x + diff.y * diff.y + diff.z * diff.z) - radius*radius;
  return beta * beta - gamma >= 0;
}

/*
Tests whether a 3D Ray hits a Triangle and returns barycentric coordinates of the intersection point
*/
__device__ bool rayIntersectsTriangle(float4 p0, float4 p1, float4 p2, float4 r, float4 o, float* u_t, float* v_t) {
  float4 edge1;
  edge1.x = p1.x - p0.x;
  edge1.y = p1.y - p0.y;
  edge1.z = p1.z - p0.z;

  float4 edge2;
  edge2.x = p2.x - p0.x;
  edge2.y = p2.y - p0.y;
  edge2.z = p2.z - p0.z;

  // Direction pvec = dir*edge2;
  float4 pvec;
  pvec.x = r.y * edge2.z - r.z * edge2.y;
  pvec.y = r.z * edge2.x - r.x * edge2.z;
  pvec.z = r.x * edge2.y - r.y * edge2.x;

  // float det = edge1%pvec;
  float det = edge1.x * pvec.x + edge1.y * pvec.y + edge1.z * pvec.z;

  if (det > -2e-6f && det < 2e-6f)
    return false;

  float inv_det = 1.0f / det;

  // 	tvec = Direction(&orig, &p0);
  float4 tvec;
  tvec.x = o.x - p0.x;
  tvec.y = o.y - p0.y;
  tvec.z = o.z - p0.z;

  // float u = (tvec%pvec) * inv_det;
  float u = (tvec.x * pvec.x + tvec.y * pvec.y + tvec.z * pvec.z) * inv_det;
  if (u < 2e-6f || u > 1.00f + 2e-6f)
    return false;

  // qvec = tvec*edge1;
  float4 qvec;
  qvec.x = tvec.y * edge1.z - tvec.z * edge1.y;
  qvec.y = tvec.z * edge1.x - tvec.x * edge1.z;
  qvec.z = tvec.x * edge1.y - tvec.y * edge1.x;

  // 	float v = (dir%qvec) * inv_det;
  float v = (r.x * qvec.x + r.y * qvec.y + r.z * qvec.z) * inv_det;
  if (v < 2e-6f || u + v > 1.0f + 2e-6f)
    return false;

  // float t = (edge2%qvec) * inv_det;
  float t = (edge2.x * qvec.x + edge2.y * qvec.y + edge2.z * qvec.z) * inv_det;
  if (t < 2e-1f)
    return false;

  *u_t = u;
  *v_t = v;
  return true;
}

__device__ void findIntersectionPoint(float4 p0, float4 p1, float4 p2, float4 n0, float4 n1, float4 n2, float u,
                                      float v, float4* sp, float4* n) {
  // interpolate point coordinates
  (*sp).x = (1 - u - v) * p0.x + u * p1.x + v * p2.x;
  (*sp).y = (1 - u - v) * p0.y + u * p1.y + v * p2.y;
  (*sp).z = (1 - u - v) * p0.z + u * p1.z + v * p2.z;

  // interpolate normal in intersection point
  (*n).x = (1 - u - v) * n0.x + u * n1.x + v * n2.x;
  (*n).y = (1 - u - v) * n0.y + u * n1.y + v * n2.y;
  (*n).z = (1 - u - v) * n0.z + u * n1.z + v * n2.z;
}

__device__ bool traceRay(float4 r, float4 o, float4* sp, float4* n, int* objID, float4* spheres, int4* faces,
                         float4* vertices, float4* normals, int* counts, const int meshcnt, bool find_ip) {
  // Für jedes Netz Schnitttests durchführen
  bool intfound = false;
  float mindist = 30E30;
  *objID        = -1;
  int absfcnt   = 0;

  for (int j = 0; j < meshcnt; j++) {
    float4 bsphere = spheres[j];
    bool indi      = rayIntersectsSphere(bsphere, r, o);
    if (indi) {
      float4 actIP;
      float4 actN;

      for (int i = 0; i < counts[j]; i++) {
        float u = 0.0f;
        float v = 0.0f;

        // Aktuelles Dreieck auslesen
        int4 t    = faces[absfcnt + i];
        float4 p0 = vertices[t.x];
        float4 p1 = vertices[t.y];
        float4 p2 = vertices[t.z];

        float4 n0 = normals[t.x];
        float4 n1 = normals[t.y];
        float4 n2 = normals[t.z];

        // Schnitttest mit Dreieck durchführen
        if (!rayIntersectsTriangle(p0, p1, p2, r, o, &u, &v)) {
          continue;
        }

        else {
          // Wenn Schnitt prüfen ob der näheste
          if (find_ip == false)
            return true;

          intfound = true;
          findIntersectionPoint(p0, p1, p2, n0, n1, n2, u, v, &actIP, &actN);
          float actdist = distanceBetweenPoints(actIP, o);

          if (actdist < mindist) {
            *objID  = j;
            mindist = actdist;
            *sp     = actIP;
            *n      = actN;
          }
        }
      }
    }

    // Offsets für Knoten und Dreiecke anpassen
    absfcnt += counts[j];
  }

  if (!intfound)
    return false;

  return true;
}

__device__ float4 shadePoint(float4 col, float4 sp, float4 n, int objID, float4 eye_r, float4* spheres, int4* faces,
                             float4* vertices, float4* normals, int* counts, const int meshcnt, const float* lightdat,
                             float* materials, float4* r, int* matidx) {
  // Material auslesen
  float mat[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (int i = 0; i < 16; i++) {
    mat[i] = materials[matidx[objID] * 16 + i];
  }

  // ambient
  col.x += lightdat_Buffer[3] * mat[3];
  col.y += lightdat_Buffer[4] * mat[4];
  col.z += lightdat_Buffer[5] * mat[5];

  n = normalizeDirection(n);

  // diffuse
  float4 light;
  light.x = lightdat_Buffer[0] - sp.x;
  light.y = lightdat_Buffer[1] - sp.y;
  light.z = lightdat_Buffer[2] - sp.z;
  light   = normalizeDirection(light);

  // scalarproduct
  float NDotL = light.x * n.x + light.y * n.y + light.z * n.z;
  if (NDotL < 0.0f)
    NDotL = 0.0f;

  eye_r = invertDirection(eye_r);
  *r    = reflect3DLA(eye_r, n);
  *r    = normalizeDirection(*r);

  if (NDotL > 0.0f) {
    // Test for Shadowrays
    float4 tmp;
    int tmpi;
    if (traceRay(light, sp, &tmp, &tmp, &tmpi, spheres, faces, vertices, normals, counts, meshcnt, false))
      return col;

    col.x += NDotL * lightdat_Buffer[6] * mat[6];
    col.y += NDotL * lightdat_Buffer[7] * mat[7];
    col.z += NDotL * lightdat_Buffer[8] * mat[8];

    // spekluar
    float RDotE = (*r).x * light.x + (*r).y * light.y + (*r).z * light.z;
    // if (RDotE < 0.0f)
    //	RDotE = 0.0f;

    if (RDotE > 0.0f) {
      col.x += lightdat_Buffer[9] * mat[9] * powf(RDotE, mat[12]);
      col.y += lightdat_Buffer[10] * mat[10] * powf(RDotE, mat[12]);
      col.z += lightdat_Buffer[11] * mat[11] * powf(RDotE, mat[12]);
    }
  }

  return col;
}

__global__ void traceLA(float4* output, const int2 resinfo, float4* spheres, int4* faces, float4* vertices,
                        float4* normals, int* counts, const int meshcnt, float* materials, int* matidx,
                        float4* objbounds, int* objmeshes, const int objcnt) {
  // Locate the actual thread in the global scope

  // TODO Aufgabe 5):
  int globalIdx      = 0;
  int globalIdy      = 0;
  int globalIdSerial = 0;

  globalIdx = blockIdx.x * blockDim.x + threadIdx.x;
  globalIdy = blockIdx.y * blockDim.y + threadIdx.y;

  if(globalIdx >= resinfo.x || globalIdy >= resinfo.y)
    return;
  
  globalIdSerial = (globalIdy * resinfo.x) + globalIdx;

  // BLOCK TO TEST INDICES

  // BLOCK TO TEST INDICES

  // Sichtstrahl für Position erstellen
  float4 eye_r = getCamRay3D(camera_Buffer, globalIdx, globalIdy, resinfo);
  float4 campos;
  campos.x = camera_Buffer[0];
  campos.y = camera_Buffer[1];
  campos.z = camera_Buffer[2];

  // Farbe des Ausgabepixels
  float4 col;
  col.x = 0.0f;
  col.y = 0.0f;
  col.z = 0.0f;

  float4 sp;  // Schnittpunkt mit der Dreiecksfläche
  float4 n;   // Interpolierte Flächennormale
  int objID = -1;

  // Schnitttest
  if (!traceRay(eye_r, campos, &sp, &n, &objID, spheres, faces, vertices, normals, counts, meshcnt, true))
    return;

  // Farbe berechnen
  float4 rdir;
  col = shadePoint(col, sp, n, objID, eye_r, spheres, faces, vertices, normals, counts, meshcnt, 0, materials, &rdir,
                   matidx);

  // Farbwert in Ausgabe
  output[globalIdSerial] = col;

  // Reflektion in der ersten Stufe berechnen
  objID = -1;
  float4 nsp;
  float4 nn;
  float4 nrdir;

  if (!traceRay(rdir, sp, &nsp, &nn, &objID, spheres, faces, vertices, normals, counts, meshcnt, true))
    return;

  col = shadePoint(col, nsp, nn, objID, rdir, spheres, faces, vertices, normals, counts, meshcnt, 0, materials, &nrdir,
                   matidx);
  output[globalIdSerial] = col;

  return;
}

cudaError_t traceWithCuda(float* result, const int blocks_x, const int blocks_y, SplitScene* scene, const int width,
                          const int height, const int facecnt, float* spheres, int* faces, float* vertices,
                          float* normals, int* loc_counts, float* materials, int* matidx, float* objbounds,
                          int* objmeshes) {
  cudaError_t cudaStatus;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Initialize Device
  cudaStatus = setupCuda();
  if (cudaStatus != cudaSuccess) {
    return cudaStatus;
  }

  // Copy Data to Device Memory
  cudaStatus = setupCudaMemory(scene, width, height, facecnt, spheres, faces, vertices, normals, loc_counts, materials,
                               matidx, objbounds, objmeshes);
  if (cudaStatus != cudaSuccess) {
    return cudaStatus;
  }

  // TODO Aufgabe 10
  // Replace this default with Aufgabe 10:

  // dim3 threads(8, 8);
  // dim3 grid(width / threads.x, height / threads.y);

  int blockSize;
  int minGridSize;

  cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, (void*) traceLA);

  cout << "Recommended blockSize is (" << blockSize << ")"
       << " and minGridSize is (" << minGridSize << ")" << endl;


  //Setzung der Blockdimensionen, sowie Grid-Dimensionen basierend auf dem Seitenverhältnis des Bildes
  float ratio = (float)width / (float)height;
  int threadsX = (int) sqrt(ratio * (float)blockSize);
  int threadsY = (int) sqrt((float)blockSize / ratio);

  dim3 threads(threadsX, threadsY);
  
  int x = ceil((float) width / (float) threadsX);
  int y = ceil((float) height / (float) threadsY);

  //Gleichmäßige Erhöhung der beiden Grid-Dimensionswerte, falls diese den Mindestwert von cudaOccupancyMaxPotentialBlockSize() unterschreiten
  if( (x * y) < minGridSize)
  {
    int rest = minGridSize - (x * y);
    int sqrt_rest = ceil(sqrt(rest));
    x += sqrt_rest;
    y += sqrt_rest;
  }
  
  dim3 grid(x, y);

  cout << "Executing Kernel with grid (" << grid.x << "," << grid.y << ")"
       << " and block (" << threads.x << "," << threads.y << ")" << endl;

  cudaEventRecord(start, 0);
  // Execute
  traceLA<<<grid, threads>>>(result_Buffer, resinfo_Buffer, spheres_Buffer, faces_Buffer, vertices_Buffer,
                             normals_Buffer, counts_Buffer, scene->getMeshCnt(), materials_Buffer, matidx_Buffer,
                             objbounds_Buffer, objmeshes_Buffer, scene->objcnt);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  // Check for any errors launching the kernel
  cudaStatus = cudaGetLastError();
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "traceLA launch failed: %s\n", cudaGetErrorString(cudaStatus));
  }

  // cudaDeviceSynchronize waits for the kernel to finish, and returns
  // any errors encountered during the launch.
  cudaStatus = cudaDeviceSynchronize();
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
  }

  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);

  std::cout << "GPU rendering required " << milliseconds / 1000.0f << "s." << std::endl;

  // Copy output from GPU to host

  // TODO Aufgabe 4)
  cudaStatus = cudaMemcpy(result, result_Buffer, width * height * sizeof(float4), cudaMemcpyDeviceToHost);
  if(cudaStatus != cudaSuccess)
  {
    fprintf(stderr, "cudaMemcpy failed! Couldn't copy \"result_Buffer\" to \"result\" from Device to Host");
    return cudaStatus;
  }

  cudaStatus = freeCudaMemory();
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "freeCudaMemory failed!");
  }

  return cudaStatus;
}