#ifndef LAGPUTRACERGLUT_H
#define LAGPUTRACERGLUT_H

#include "GA_GPU_Tracer_glut.h"

class LA_GPU_Tracer : public GA_GPU_Tracer
{
public:
	LA_GPU_Tracer(int width, int height) : GA_GPU_Tracer(width, height)
	{
		blockSize = 8;
		this->height = height;
		this->width = width;
		this->type=1;
	}

LA_GPU_Tracer() : GA_GPU_Tracer(width, height)  
	{
		blockSize = 8;
		this->height = height;
		this->width = width;
		this->type=1;
    }
};

#endif


