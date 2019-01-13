#include "material.h"

Material::Material()
{
	ambient = Vector(0.0f);
	diffuse = Vector(0.0f);
	specular = Vector(0.0f);
	mirror = Vector(0.0f);
	transparency = Vector(0.0f);

	roughness = 0;
	refraction_index = 1;
	phong_exponent = 0;

	brdf_id = -1;
}