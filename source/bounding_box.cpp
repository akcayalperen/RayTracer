#include "bounding_box.h"
#include "scene.h"
#include <algorithm>

BoundingBox::BoundingBox()
{
	min = Vector();
	max = Vector();
}

BoundingBox::BoundingBox(Vector min, Vector max)
{
	this->min = min;
	this->max = max;
}

void BoundingBox::transformBBox(Matrix matrix, Vector motion_vector)
{
	Vector corners[8] = {
		Vector(min.x, min.y, min.z),
		Vector(max.x, min.y, min.z),
		Vector(min.x, max.y, min.z),
		Vector(min.x, min.y, max.z),
		Vector(max.x, max.y, min.z),
		Vector(max.x, min.y, max.z),
		Vector(min.x, max.y, max.z),
		Vector(max.x, max.y, max.z),
	};

	Vector results[8];

	for(int i = 0; i < 8; i++)
	{
		results[i] = matrix.transform_point(corners[i]);
	}

	for (int i = 0; i < 8; i++)
	{
		min.x = std::min(results[i].x, min.x);
		min.y = std::min(results[i].y, min.y);
		min.z = std::min(results[i].z, min.z);

		max.x = std::max(results[i].x, max.x);
		max.y = std::max(results[i].y, max.y);
		max.z = std::max(results[i].z, max.z);
	}

	if (motion_vector.x + motion_vector.y + motion_vector.z > 0.0f)
	{
		Matrix motion_matrix = Matrix::TranslationMatrix(motion_vector);
		
		motion_matrix = motion_matrix * matrix;

		for (int i = 0; i < 8; i++)
		{
			results[i] = motion_matrix.transform_point(corners[i]);
		}

		for (int i = 0; i < 8; i++)
		{
			min.x = std::min(results[i].x, min.x);
			min.y = std::min(results[i].y, min.y);
			min.z = std::min(results[i].z, min.z);

			max.x = std::max(results[i].x, max.x);
			max.y = std::max(results[i].y, max.y);
			max.z = std::max(results[i].z, max.z);
		}

	}
}

bool BoundingBox::intersect(const Ray &ray)
{
	float tmin = (min.x - ray.origin.x) / ray.direction.x;
	float tmax = (max.x - ray.origin.x) / ray.direction.x;

	if (tmin > tmax) std::swap(tmin, tmax);

	float tymin = (min.y - ray.origin.y) / ray.direction.y;
	float tymax = (max.y - ray.origin.y) / ray.direction.y;

	if (tymin > tymax) std::swap(tymin, tymax);

	if ((tmin > tymax) || (tymin > tmax))
		return false;

	if (tymin > tmin)
		tmin = tymin;

	if (tymax < tmax)
		tmax = tymax;

	float tzmin = (min.z - ray.origin.z) / ray.direction.z;
	float tzmax = (max.z - ray.origin.z) / ray.direction.z;

	if (tzmin > tzmax) std::swap(tzmin, tzmax);

	if ((tmin > tzmax) || (tzmin > tmax))
		return false;

	if (tzmin > tmin)
		tmin = tzmin;

	if (tzmax < tmax)
		tmax = tzmax;

	return true;
}