#pragma once
#ifndef __OBJECT__H_
#define __OBJECT__H_

#include "ray.h"
#include "rtMatrix.h"
#include "bounding_box.h"
#include "intersection.h"
#include <vector>


class Object
{
public:
	enum class ShadingMode
	{
		Flat,
		Smooth
	};

	Object()
	{

	}

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) = 0;

	ShadingMode shading_mode;
	BoundingBox bbox;
	Vector motion_vector;

	Matrix base_transformation;
	Matrix inverse_transformation;
	Matrix transpose_transformation;

};

#endif