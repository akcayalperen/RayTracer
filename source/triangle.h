#pragma once
#ifndef __TRIANGLE__H_
#define __TRIANGLE__H_

#include "object.h"

class Triangle : public Object
{
public:
	Triangle();
	Triangle(int v0, int v1, int v2, int vertex_offset);

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;
	void calculate_dp();

	int v0_id;
	int v1_id;
	int v2_id;

	int vertex_offset;
	int texture_offset;
	Vector dpdu;
	Vector dpdv;
};

#endif