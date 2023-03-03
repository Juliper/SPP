#ifndef VIRTUALTRACERGLUT_H
#define VIRTUALTRACERGLUT_H

#include "Color.h"
#include "Material.h"
#include "Light.h"
#include "SplitScene.h"
#include "Constants.h"

class virtualTracer
{

public:

	virtualTracer(int width, int height);
    ~virtualTracer();

	void renderImage();

	float* checkImage;

	bool rendered;

	int height;
	int width;

	int facecnt;
	int vtxcnt;
	int objcnt;

	float * spheres;
	int * faces;
	int * faces1;
	int * matidx;

	float * vertices;
	float * vertices1;

	float * vertices_norm;
	float * vertices_norm1;

	int * loc_counts;

	float * materials;
	float * materials1;

	float * objbounds;
	int * objmeshes;

	SplitScene* scene;
};

#endif