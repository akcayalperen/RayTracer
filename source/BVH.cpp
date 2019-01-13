#include "BVH.h"
#include "scene.h"

BoundingBox recalculate_bbox(const BoundingBox &b1, const BoundingBox &b2)
{
	float minx = std::fmin(b1.min.x, b2.min.x);
	float miny = std::fmin(b1.min.y, b2.min.y);
	float minz = std::fmin(b1.min.z, b2.min.z);

	float maxx = std::fmax(b1.max.x, b2.max.x);
	float maxy = std::fmax(b1.max.y, b2.max.y);
	float maxz = std::fmax(b1.max.z, b2.max.z);

	return BoundingBox(Vector(minx, miny, minz), Vector(maxx, maxy, maxz));
}

BVH::BVH(std::vector<Object*> &objects, SplitRange range, int axis)
{
	bbox = objects[range.start_index]->bbox;
	for (int i = range.start_index; i < range.end_index; i++)
	{
		bbox = recalculate_bbox(bbox, objects[i]->bbox);
	}

	if (objects.size() == 1)
	{
		left = objects[0];
		right = objects[0];
	}
	else if (objects.size() == 2)
	{
		left = objects[0];
		right = objects[1];
	}
	else
	{
		Vector split = (bbox.max + bbox.min) / 2.0f;

		int result = find_split_range(objects, split.coordinate[axis], axis, range);

		SplitRange left_range = { range.start_index, result };
		SplitRange right_range = { result, range.end_index };

		left = create_branch(objects, left_range, (axis + 1) % 3);
		right = create_branch(objects, right_range, (axis + 1) % 3);
	}
}

int BVH::find_split_range(std::vector<Object*> &objects, float split, int axis, SplitRange range)
{
	int index = range.start_index;

	for (int i = range.start_index; i < range.end_index; i++)
	{
		float center = (objects[i]->bbox.min.coordinate[axis] + objects[i]->bbox.max.coordinate[axis]) / 2.0f;

		if (center < split)
		{
			std::swap(objects[i], objects[index]);
			index++;
		}
	}

	if (index == range.start_index || index == range.end_index - 1)
	{
		index = (range.start_index + range.end_index) / 2;
	}

	return index;
}

Object* BVH::create_branch(std::vector<Object*> &objects, SplitRange range, int axis)
{
	if (range.end_index - range.start_index == 1)
	{
		return objects[range.start_index];
	}

	return new BVH(objects, range, axis);
}

bool BVH::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	if (!bbox.intersect(ray))
	{
		return false;
	}

	bool left_intersect = left->intersect(ray, intersection, culling, texture);
	bool right_intersect = right->intersect(ray, intersection, culling, texture);

	return (left_intersect || right_intersect);
}