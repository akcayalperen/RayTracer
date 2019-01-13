#include "BRDF.h"
#include "scene.h"
#include <algorithm>
#include <cmath>
#include <iostream>

BRDF::BRDF()
{
	exponent = 0;
	normalizer = 0;
}

//Original Phong
OriginalPhongBRDF::OriginalPhongBRDF()
{
	exponent = 0;
	normalizer = 0;
}

Vector OriginalPhongBRDF::shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector)
{
	Vector light_direction = light_vector.normalized();
	Vector surface_normal = intersection.surface_normal;

	float cosinetheta = std::max(0.0f, light_direction.dot(surface_normal));

	if (cosinetheta < 90 * PI / 180)
	{
		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector w_r = (-light_direction + surface_normal * 2 * surface_normal.dot(light_direction)).normalized();
		float cosinealpha = std::max(0.0f, w_r.dot(-normalized_direction));

		return (rt_scene->materials[intersection.material_id].diffuse + rt_scene->materials[intersection.material_id].specular *
			pow(cosinealpha, exponent) / cosinetheta) * cosinetheta * light->calculate_intensity(light_vector, intersection.pw);
	}

	return Vector(0.0f);
}

//Modified Phong
ModifiedPhongBRDF::ModifiedPhongBRDF()
{
	exponent = 0;
	normalizer = 0;
}

Vector ModifiedPhongBRDF::shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector)
{
	Vector light_direction = light_vector.normalized();
	Vector surface_normal = intersection.surface_normal;

	float cosinetheta = std::max(0.0f, light_direction.dot(surface_normal));

	if (cosinetheta < 90 * PI / 180)
	{
		bool normalized = normalizer != 0;

		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector w_r = (-light_direction + surface_normal * 2 * surface_normal.dot(light_direction)).normalized();
		float cosinealpha = std::max(0.0f, w_r.dot(-normalized_direction));

		Vector diffuse = rt_scene->materials[intersection.material_id].diffuse;
		Vector specular = rt_scene->materials[intersection.material_id].specular;

		if (normalized)
		{
			diffuse = diffuse / PI;
			specular = specular * ((exponent + normalizer) / (normalizer * PI));
		}

		return (diffuse + specular *
			pow(cosinealpha, exponent)) * cosinetheta * light->calculate_intensity(light_vector, intersection.pw);
	}

	return Vector(0.0f);
}

//Original Blinn-Phong
OriginalBlinnPhongBRDF::OriginalBlinnPhongBRDF()
{
	exponent = 0;
	normalizer = 0;
}

Vector OriginalBlinnPhongBRDF::shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector)
{
	Vector light_direction = light_vector.normalized();
	Vector surface_normal = intersection.surface_normal;

	float cosinetheta = std::max(0.0f, light_direction.dot(surface_normal));

	if (cosinetheta < 90 * PI / 180)
	{
		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector half_vector = (light_direction - normalized_direction).normalized();
		float cosinealpha = std::max(0.0f, half_vector.dot(surface_normal));

		return (rt_scene->materials[intersection.material_id].diffuse + rt_scene->materials[intersection.material_id].specular *
			pow(cosinealpha, exponent) / cosinetheta) * cosinetheta * light->calculate_intensity(light_vector, intersection.pw);
	}

	return Vector(0.0f);
}

//Modified Blinn-Phong
ModifiedBlinnPhongBRDF::ModifiedBlinnPhongBRDF()
{
	exponent = 0;
	normalizer = 0;
}

Vector ModifiedBlinnPhongBRDF::shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector)
{
	Vector light_direction = light_vector.normalized();
	Vector surface_normal = intersection.surface_normal;

	float cosinetheta = std::max(0.0f, light_direction.dot(surface_normal));

	if (cosinetheta < 90 * PI / 180)
	{
		bool normalized = normalizer != 0;

		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector half_vector = (light_direction - normalized_direction).normalized();
		float cosinealpha = std::max(0.0f, half_vector.dot(surface_normal));

		Vector diffuse = rt_scene->materials[intersection.material_id].diffuse;
		Vector specular = rt_scene->materials[intersection.material_id].specular;

		if (normalized)
		{
			diffuse = diffuse / PI;
			specular = specular * ((exponent + normalizer) / (normalizer * PI));
		}

		return (diffuse + specular *
			pow(cosinealpha, exponent)) * cosinetheta * light->calculate_intensity(light_vector, intersection.pw);
	}

	return Vector(0.0f);
}

//Torrance Sparrow
TorranceSparrowBRDF::TorranceSparrowBRDF()
{
	exponent = 0;
	normalizer = 0;
	refractive_index = 0;
}

Vector TorranceSparrowBRDF::shade(const Ray& ray, const Intersection& intersection, Light *light, Vector light_vector)
{
	Vector light_direction = light_vector.normalized();
	Vector surface_normal = intersection.surface_normal;

	float cosinetheta = std::max(0.0f, light_direction.dot(surface_normal));

	if (cosinetheta < 90 * PI / 180)
	{
		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector half_vector = (light_direction - normalized_direction).normalized();
		float cosinealpha = std::max(0.0f, half_vector.dot(surface_normal));

		float probability = ((exponent + 2) / (2 * PI)) * pow(cosinealpha, exponent);

		float wo_dot_wh = (-normalized_direction).dot(half_vector);

		float geometry = std::min(1.0f, 
			std::min(2 * surface_normal.dot(half_vector) * surface_normal.dot(-normalized_direction) / wo_dot_wh,
			2 * surface_normal.dot(half_vector) * surface_normal.dot(light_direction) / wo_dot_wh));

		float ro = ((refractive_index - 1) * (refractive_index - 1)) / ((refractive_index + 1) * (refractive_index + 1));

		float cosinebeta = std::max(0.0f, half_vector.dot(-normalized_direction));

		float fresnel = ro + (1 - ro) * pow(1 - cosinebeta, 5);

		float cosinephi = surface_normal.dot(-normalized_direction);

		Vector diffuse = rt_scene->materials[intersection.material_id].diffuse / PI;
		Vector specular = rt_scene->materials[intersection.material_id].specular;

		Vector result = (diffuse + 
			specular * (probability * fresnel * geometry) / (4 * cosinetheta * cosinephi)) *
			cosinetheta * light->calculate_intensity(light_vector, intersection.pw);

		return result;
	}

	return Vector(0.0f);
}