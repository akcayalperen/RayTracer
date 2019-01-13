#pragma once
#ifndef __MESH__H_
#define __MESH__H_

#include "BVH.h"
#include "object.h"
#include "triangle.h"
#include "light.h"
#include <vector>

class Mesh : public Object
{
public:

	Mesh();

	void construct_bvh();
	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;

	int material_id;
	int texture_id;
	std::vector<Object*> faces;
	
	BVH *bvh;
};

class MeshInstance : public Object
{
public:

	MeshInstance();

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;

	int material_id;
	int texture_id;
	int base_mesh_id;
};

class LightMesh : public Mesh, public Light
{
public:
	LightMesh();
	void calculate_area();

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;
	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	std::vector<float> areas;
	Vector radiance;
	float total_area;
};

#endif