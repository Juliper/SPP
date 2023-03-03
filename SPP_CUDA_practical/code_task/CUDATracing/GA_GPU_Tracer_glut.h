#ifndef GAGPUTRACERGLUT_H
#define GAGPUTRACERGLUT_H

#include "virtualTracer_glut.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "timing.h"

class GA_GPU_Tracer : public virtualTracer//, public SDKSample
{
public:
	GA_GPU_Tracer(int width, int height) : virtualTracer(width, height)    
	{
		blockSize = 8;
		this->height = height;
		this->width = width;
		this->type=0;
		this->checkImage = new float [width*height*3];
    }

	//int setupMatrixTranspose(); //
    int setupCL(); //
    int runCLKernels(); //
    int initialize(); //
    int setup(); //
	int run();//
	int cleanup();//
	int verifyResults(); //
	void printStats(); //
	void renderImage();

	int cleanupCL(void); //

	void renderWithOpenGL();

	/*cl_mem spheres_Buffer;
	cl_mem objbounds_Buffer;
	cl_mem objmeshes_Buffer;
	cl_mem faces_Buffer;
	cl_mem vertices_Buffer;
	cl_mem loc_counts_Buffer;
	cl_mem counts_Buffer;
	cl_mem materials_Buffer;
	cl_mem matidx_Buffer;
	cl_mem result_Buffer;
	cl_mem result_AnaVals;*/

protected:
	int type;
	unsigned int blockSize;
	int width;
	int height;

	/*cl_platform_id		platform;
	cl_context          context;
	cl_device_id        *devices;
	cl_command_queue    commandQueue;

	cl_program program;

	// This program uses only one kernel and this serves as a handle to it
	cl_kernel  kernel;*/

private:
	std::string convertToString(const char *filename);
};

#endif