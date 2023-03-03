#ifndef SPLITSCENE_H
#define SPLITSCENE_H

#include <vector>
#include "MeshSphere.h"
#include "DirectionPoint.h"
#include "Camera.h"
#include "Light.h"
#include "Material.h"
//#include <fstream>
//#include <string>
//#include <cstring>
//#include <iostream>
#include "Log.h"

using namespace std;

class Mesh;
class Mesh2;
class Sphere;

class SplitScene
{
public:
	static SplitScene* Instance(); 

	Light light;
	Point opengllight;
	Camera cam;

	SplitScene();
	void addMesh(Mesh m);
	void addMaterial(Material mat);
	void configureScene();
	void clearScene();

	int getMaterialCount();
	int getMeshCnt();
	int getOldMeshCnt();
	int getSplitFac();
	void setSplitFac(int fac);

	std::vector<Mesh2> meshes;
	std::vector<Mesh2> obj;

	std::vector<Mesh> oldmeshes;
	std::vector<Point> vtx;
	std::vector<Point> openglvtx;
	std::vector<Material> mats;

	int vtxcnt;
	int openglvtxcnt;

	bool force_cpu_fallback;

	// Daten für die Objekte der Szene
	// Die gebraucht werden, wenn mit den großen Bounds gearbeitet wird
	int objcnt;
	std::vector<Sphere> objbounds;
	std::vector<int> objmeshes;

	void rotateScene(float angle, float x, float y, float z);
	void resetScene();
	void TranslateSceneRTM();
	void TranslateSceneMEM();
	void exportTriangle(int i);

private:
	static SplitScene* _instance;

	int meshcnt;
	int oldmeshcnt;
	int splitfac;
	int matcount;
	ofstream rtm_stream;
	ofstream mem_stream;
	ofstream tria_stream;

	void splitMesh(Mesh2 m);
	void addConfiguredMesh(Mesh2 m);
	void WriteRTMLine(const char* logline,...);
	void WriteMEMLine(const char* logline,...);
	void WriteTriaLine(const char* logline,...);
};

#endif