#ifndef RTMLOADER_H
#define RTMLOADER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MeshSphere.h"
#include "DirectionPoint.h"
#include "Transformation.h"
#include "Face.h"

class rtmLoader {
public:
	rtmLoader();
	Mesh loadRtm(const std::string &filename, Transformation* trans);	
};
#endif