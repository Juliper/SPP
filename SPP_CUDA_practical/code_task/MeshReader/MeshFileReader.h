#ifndef MESHFILEREADER_H
#define MESHFILEREADER_H

#include "objLoader.h"
#include "rtmLoader.h"
#include "MeshSphere.h"
#include "DirectionPoint.h"
#include "Transformation.h"

class MeshFileReader
{
public:
	MeshFileReader();
	Mesh readMeshFromFile(char *filename, Transformation* trans, int filetype);	 

private:
	Mesh readObjFile(char *filename, Transformation* trans);
	Mesh readRtmFile(char *filename, Transformation* trans);
};
#endif