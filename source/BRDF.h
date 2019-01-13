#pragma once
#ifndef __BRDF__H_
#define __BRDF__H_
#include "rtVector.h"
#include "light.h"
#include "ray.h"
#include "intersection.h"

class BRDF
{
public:
	BRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) = 0;

	float exponent;
	int normalizer;
};

class OriginalPhongBRDF : public BRDF
{
public:
	OriginalPhongBRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) override;
};

class ModifiedPhongBRDF : public BRDF
{
public:
	ModifiedPhongBRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) override;
};

class OriginalBlinnPhongBRDF : public BRDF
{
public:
	OriginalBlinnPhongBRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) override;
};

class ModifiedBlinnPhongBRDF : public BRDF
{
public:
	ModifiedBlinnPhongBRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) override;
};

class TorranceSparrowBRDF : public BRDF
{
public:
	TorranceSparrowBRDF();
	virtual Vector shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector) override;

	float refractive_index;
};

#endif // !__BRDF__H_
