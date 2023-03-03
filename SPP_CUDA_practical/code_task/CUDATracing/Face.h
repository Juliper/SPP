#ifndef FACE_H
#define FACE_H
#include "DirectionPoint.h"

class Face {
public:
	int idx;
	int p0;
	int p1;
	int p2;

	double u_t;
	double v_t;

	Face (int p0, int p1, int p2);
	Face() {};
};

#endif