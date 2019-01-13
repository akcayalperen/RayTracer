#pragma once
#ifndef __RAY__H_
#define __RAY__H_

#include "rtVector.h"

class Ray
{
public:
	Ray(Vector origin, Vector direction);

	Vector origin;
	Vector direction;

	float delta_time;
	bool inside;

	bool light_hit;

	float target_pixel_x;
	float target_pixel_y;
};

#endif