#ifndef KERNEL_CUH
#define KERNEL_CUH

#include "SplitScene.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

cudaError_t addWithCuda(int* c, const int* a, const int* b, unsigned int size);
cudaError_t traceWithCuda(float* result, const int blocks_x, const int blocks_y, SplitScene* scene, const int width,
                          const int height, const int facecnt, float* spheres, int* faces, float* vertices,
                          float* normals, int* loc_counts, float* materials, int* matidx, float* objbounds,
                          int* objmeshes);

__global__ void traceLA(float4* output, const int2 resinfo,
                        /*const	  float* camera,*/
                        float4* spheres, int4* faces, float4* vertices, float4* normals, int* counts, const int meshcnt,
                        /*const float* lightdat,*/
                        float* materials, int* matidx, float4* objbounds, int* objmeshes, const int objcnt);

#endif