#pragma once
#ifndef __INTERSECTION__H_
#define __INTERSECTION__H_
#include "rtVector.h"
#include "texture.h"

class Intersection
{
public:
	Intersection();

	float t;
	Vector surface_normal;
	int material_id;

	float u;
	float v;
	Vector texture_color;
	Texture::DecalMode decal_mode;

	Vector dpdu;
	Vector dpdv;

	Vector radiance;

	float pw;
};

#endif // !__INTERSECTION__H_
