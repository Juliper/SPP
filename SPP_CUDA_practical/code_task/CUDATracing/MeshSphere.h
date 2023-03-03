#ifndef MESHSPHERE_H
#define MESHSPHERE_H

#include <iostream>
#include <vector>
#include <math.h>
#include "Face.h"
#include "DirectionPoint.h"
#include "Material.h"
#include "Pointpair.h"
#include "Triangle.h"
#include "Ray3D.h"
#include "SplitScene.h"

using namespace std;

class Mesh;
class Mesh2;
class Sphere;
class BoundSphereCreator;

class Sphere {
public:
	float s2;
	float s3;
	float s4;
	float s5;
	float s6;

	Point MP;
	float radius;

	Sphere();
	Sphere(Point p, float radius);
	Sphere(float s2, float s3, float s4, float s5);
	Sphere(Mesh* m);
	Sphere(Mesh2* m, std::vector<Point>  points);

	Point* giveMP();

	bool doesRayIntersect(Ray* r);
	bool doesRayIntersect2(Ray* r);
	bool doesRayIntersect(Ray3D* r);
	Point intersectSphere(Ray* r);

private:
	float distOfPoints(Point* p1, Point* p2);
};

class Mesh {
public:
	int getObjectCount();	

	int objectcount;
	int facecount;
	int vtxcount;

	vector<Face> faces;
	vector<Point> points;
	vector<Triangle> objects;
	vector<Direction> facenormals;

	Sphere bound;

	void setFaceCount(int c);
	void setVtxCount(int c);

	void addFace(int idx, Face f);
	void addPoint(int idx, Point p);
	void addFace(Face f);
	void addPoint(Point p);
	void addObject(Triangle t);

	void interpolateNormals();

	Mesh();
	void configureMesh();

	int mat;
};

class Mesh2 {
public:
	int getObjectCount();	
	int facecount;
	int vtxcount;
	int objectcount;

	vector<Face> faces;
	vector<int> points;
	vector<Triangle> objects;

	Sphere bound;

	void addFace(Face f);
	void addPoints(Face * f);
	void addObject(Triangle t);

	Mesh2();
	void configureMesh(vector<Point> points);

	int mat;

private:
	bool pointIsNew(int p);
};

class BoundSphereCreator {
public:
	BoundSphereCreator();
	float* createBound(Mesh2* m, std::vector<Point>  vtx);
};

#endif
