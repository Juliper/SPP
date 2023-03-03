#include "kernel.cuh"

/////////////// HOST FUNCTIONS ////////////////////////

cudaError_t setupCuda() {
  cudaError_t cudaStatus;
  cudaStatus = cudaSetDevice(0);
  if (cudaStatus != cudaSuccess) {
    fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
  }
  return cudaStatus;
}

/////////////// DEVICE FUNCTIONS ////////////////////////
__device__ float4 normalizeDirection(float4 d) {
  float4 nd;
  float abs = sqrtf(d.x * d.x + d.y * d.y + d.z * d.z);
  nd.x      = d.x / abs;
  nd.y      = d.y / abs;
  nd.z      = d.z / abs;
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
