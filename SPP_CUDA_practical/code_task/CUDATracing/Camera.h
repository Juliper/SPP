#ifndef CAMERA_H
#define CAMERA_H

#include "DirectionPoint.h"
#include "Ray.h"
#include "Ray3D.h"

class Point;
class Ray;
class Ray3D;

class Camera {
public:
	float cameraFOVwidth;
	float cameraFOVheight;
	int imageWidth;
	int imageHeight;
	Point* pos;

	float pixelWidth;
	float pixelHeight;

	Camera(Point* pos, float cameraFOVwidth, float cameraFOVheight, int imageWidth, int imageHeight);
	Camera();
	Ray getCamRay(int x, int y);
	Ray3D getCamRay3D(int x, int y);
};

#endif