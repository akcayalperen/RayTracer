#pragma once
#ifndef __MATERIAL__H_
#define __MATERIAL__H_

#include "rtVector.h"

class Material
{
public:
	Material();

	Vector ambient;
	Vector diffuse;
	Vector specular;
	Vector mirror;
	Vector transparency;

	float roughness;
	float phong_exponent;
	float refraction_index;

	int brdf_id;
};

#endif