#pragma once
#ifndef __BVH__H_
#define __BVH__H_

#include "object.h"
#include <vector>

class BVH : public Object
{
public:
	struct SplitRange
	{
		int start_index;
		int end_index;
	};

	BVH(std::vector<Object*> &objects, SplitRange range, int axis = 0);

	int find_split_range(std::vector<Object*> &objects, float split, int axis, SplitRange range);
	Object* create_branch(std::vector<Object*> &objects, SplitRange range, int axis = 0);

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;

	Object* left;
	Object* right;
};

#endif