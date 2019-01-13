#pragma once
#ifndef __LIGHT__H_
#define __LIGHT__H_

#include "rtVector.h"
#include "ray.h"
#include "intersection.h"

class Light
{
public:
	Vector position;

	virtual Vector calculate_intensity(Vector& light_vector, float cosine) = 0;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) = 0;
};


class PointLight : public Light
{
public:
	PointLight();

	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	Vector intensity;

};

class AreaLight : public Light
{
public:
	AreaLight();

	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	Vector intensity;

	Vector edge_vector_1;
	Vector edge_vector_2;
	Vector normal;
};

class DirectionalLight : public Light
{
public:
	DirectionalLight();

	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	Vector direction;
	Vector radiance;
};

class SpotLight : public Light
{
public:
	SpotLight();

	virtual Vector calculate_intensity(Vector& light_vector, float cosine) override;
	virtual Vector get_position(const Ray& ray, Intersection& intersection) override;

	Vector direction;
	Vector intensity;

	float coverage_angle;
	float falloff_angle;
};

#endif