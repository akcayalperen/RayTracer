#include "intersection.h"

#define FAR_FAR_AWAY 999999.0f

Intersection::Intersection()
{
	surface_normal = Vector(0.0f);
	texture_color = Vector(0.0f);
	decal_mode = Texture::DecalMode::None;
	t = FAR_FAR_AWAY;
	radiance = Vector(0.0f);
}