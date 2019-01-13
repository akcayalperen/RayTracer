#pragma once
#ifndef __BOUNDING_BOX__H_
#define __BOUNDING_BOX__H_

#include "ray.h"
#include "rtMatrix.h"

class BoundingBox
{
public:
	BoundingBox();
	BoundingBox(Vector min, Vector max);

	void transformBBox(Matrix matrix, Vector motion_vector);
	bool intersect(const Ray& ray);

	Vector min;
	Vector max;
};

#endif