#ifndef MATERIALBUILDER_H
#define MATERIALBUILDER_H
#include "Material.h"

class MaterialBuilder
{
public:
	MaterialBuilder();
	static MaterialBuilder* Instance(); 
	Material getMaterial(int mat);

private:
	static MaterialBuilder* _instance;
};

#endif