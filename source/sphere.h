#pragma once
#ifndef __SPHERE__H_
#define __SPHERE__H_

#include "object.h"
#include "light.h"

class Sphere : public Object
{
public:
	Sphere();

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;


	int material_id;
	int texture_id;
	int center_vertex_id;
	float radius;
};

class LightSphere : public Sphere, public Light
{
public:
	LightSphere();

	virtual bool intersect(const Ray& ray, Intersection& intersection, bool culling = true, bool texture = false) override;
	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	Vector radiance;
};

#endif