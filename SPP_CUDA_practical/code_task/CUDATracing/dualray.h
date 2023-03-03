#ifndef DUALRAY_H
#define DUALRAY_H
#include "DirectionPoint.h"

class DualRay {
public:
	float ray18;
	float ray20;
	float ray22;
	float ray23;
	float ray25;
	float ray26;
	Point origin;

	DualRay(Point* p0, Point* p1);
};

#endif