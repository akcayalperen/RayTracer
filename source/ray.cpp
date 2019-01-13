#include "ray.h"

Ray::Ray(Vector origin, Vector direction)
{
	this->origin = origin;
	this->direction = direction;

	delta_time = 0.0f;
	inside = false;
	light_hit = false;
}