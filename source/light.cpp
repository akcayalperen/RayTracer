#include "light.h"
#include "scene.h"

PointLight::PointLight()
{
	position = Vector(0.0f);
	intensity = Vector(0.0f);
}

Vector PointLight::get_position(const Ray& ray, Intersection& intersection)
{
	return position;
}

Vector PointLight::calculate_intensity(Vector& light_vector, float cosine)
{
	return intensity / light_vector.length_squared();
}

AreaLight::AreaLight()
{
	position = Vector(0.0f);
	intensity = Vector(0.0f);
	edge_vector_1 = Vector(0.0f);
	edge_vector_2 = Vector(0.0f);
	normal = Vector(0.0f);
}

Vector AreaLight::get_position(const Ray& ray, Intersection& intersection)
{
	double r1 = rt_scene->dist(rt_scene->rng);
	double r2 = rt_scene->dist(rt_scene->rng);

	return position + edge_vector_1 * r1 + edge_vector_2 * r2;
}

Vector AreaLight::calculate_intensity(Vector& light_vector, float cosine)
{
	return (intensity * ((-light_vector).normalized()).dot(normal)) / light_vector.length_squared();
}

DirectionalLight::DirectionalLight()
{
	position = Vector(0.0f);
	direction = Vector(0.0f);
	radiance = Vector(0.0f);
}

Vector DirectionalLight::get_position(const Ray& ray, Intersection& intersection)
{
	return position;
}

Vector DirectionalLight::calculate_intensity(Vector& light_vector, float cosine)
{
	return radiance;
}

SpotLight::SpotLight()
{
	position = Vector(0.0f);
	direction = Vector(0.0f);
	intensity = Vector(0.0f);

	falloff_angle = 0.0f;
	coverage_angle = 0.0f;
}

Vector SpotLight::get_position(const Ray& ray, Intersection& intersection)
{
	return position;
}

Vector SpotLight::calculate_intensity(Vector& light_vector, float cosine)
{
	Vector normalized_light_vector = (-light_vector).normalized();

	float costheta = normalized_light_vector.dot(direction);
	float cosalpha = cos(PI * (falloff_angle / 2) / 180);
	float cosbeta = cos(PI * (coverage_angle / 2) / 180);

	if (costheta > cosalpha)
	{
		return intensity / light_vector.length_squared();
	}
	else if (costheta < cosbeta)
	{
		return Vector(0.0f);
	}
	else
	{
		float c = std::pow((costheta - cosbeta) / (cosalpha - cosbeta), 4);

		return (intensity * c) / light_vector.length_squared();
	}
}

