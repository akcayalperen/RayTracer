#include "sphere.h"
#include "ray.h"
#include "scene.h"
#include <iostream>
#include <algorithm>

#define INTERSECTION_TEST_EPSILON 0

Sphere::Sphere()
{
	motion_vector = Vector(0.0f);
	base_transformation = Matrix::Identity();
	texture_id = -1;
}

bool Sphere::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	Matrix inverse = inverse_transformation;
	Matrix transpose = transpose_transformation;

	if (motion_vector.x + motion_vector.y + motion_vector.z > 0.0f)
	{
		Matrix motion_matrix = Matrix::TranslationMatrix(motion_vector * ray.delta_time);

		Matrix transformation = motion_matrix * base_transformation;
		inverse = transformation.inverse();
		transpose = inverse.transpose();
	}

	Ray transformed_ray = Ray(inverse.transform_point(ray.origin), inverse.transform_vector(ray.direction));
	transformed_ray.inside = ray.inside;

	Vector center = rt_scene->vertex_data[center_vertex_id];

	float root1;
	float root2;

	Vector direction = transformed_ray.direction;
	Vector distance_vector = transformed_ray.origin - center;

	float A = direction.dot(direction);
	float B = direction.dot(distance_vector) * 2;
	float C = distance_vector.dot(distance_vector) - radius * radius;

	float discriminant = B * B - 4 * A * C;
	if (discriminant < INTERSECTION_TEST_EPSILON)
	{
		return 0.0f;
	}
	else
	{
		discriminant = sqrt(discriminant);
		root1 = (-B + discriminant) / (2 * A);
		root2 = (-B - discriminant) / (2 * A);
	}

	float small = root1 < root2 ? root1 : root2;
	float big = root1 > root2 ? root1 : root2;

	if (small < INTERSECTION_TEST_EPSILON)
	{
		if (big < INTERSECTION_TEST_EPSILON)
		{
			return false;
		}
		else
		{
			small = big;
		}
	}

	if (intersection.t > small)
	{

		intersection.t = small;
		intersection.material_id = material_id;
		intersection.surface_normal = (transformed_ray.origin + transformed_ray.direction * small - center).normalized();

		if (texture_id != -1)
		{
			Vector intersection_point = transformed_ray.origin + transformed_ray.direction * small - center;

			float teta = acos(intersection_point.y / radius);
			float fi = atan2(intersection_point.z, intersection_point.x);

			intersection.u = (-fi + PI) / (2 * PI);
			intersection.v = teta / PI;


			if (rt_scene->textures[texture_id].bump_map)
			{
				Vector dpdu = Vector(2 * PI * intersection_point.z, 0.0f, -2 * PI * intersection_point.x);
				Vector dpdv = Vector(PI * intersection_point.y * cos(fi), -PI * radius * sin(teta), PI * intersection_point.y * sin(fi));

				if (rt_scene->textures[texture_id].texture_type == Texture::TextureType::Perlin)
				{
					Vector gradient = rt_scene->textures[texture_id].get_perlin_bump(intersection_point.x, intersection_point.y, intersection_point.z);

					Vector gradient_parallel = intersection.surface_normal * (gradient.dot(intersection.surface_normal));
					Vector gradient_perpendicular = gradient - gradient_parallel;

					intersection.surface_normal = (intersection.surface_normal - gradient_perpendicular).normalized();
				}
				else
				{

					Vector gradient = rt_scene->textures[texture_id].get_bump(intersection.u, intersection.v);

					Vector qu = dpdu + gradient.x * intersection.surface_normal;
					Vector qv = dpdv + gradient.y * intersection.surface_normal;

					intersection.surface_normal = (qv.cross(qu)).normalized();
				}
			}

			//intersection_point = ray.origin + ray.direction * small;
			if (rt_scene->textures[texture_id].texture_type == Texture::TextureType::Perlin)
			{
				intersection.texture_color = rt_scene->textures[texture_id].get_perlin_color(intersection_point.x, intersection_point.y, intersection_point.z);
			}
			else
			{
				intersection.texture_color = rt_scene->textures[texture_id].get_color(intersection.u, intersection.v) / rt_scene->textures[texture_id].normalizer;
			}
			intersection.decal_mode = rt_scene->textures[texture_id].decal_mode;
		}
		else
		{
			intersection.decal_mode = Texture::DecalMode::None;
		}

		intersection.surface_normal = transpose.transform_vector(intersection.surface_normal).normalized();
		intersection.radiance = Vector(0.0f);

		return true;
	}

	return false;

}


//Light Sphere
LightSphere::LightSphere()
{
	motion_vector = Vector(0.0f);
	base_transformation = Matrix::Identity();

	radiance = Vector(0.0f);
}

bool LightSphere::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	bool result = false;

	if (Sphere::intersect(ray, intersection, culling, texture))
	{
		result = true;
		intersection.radiance = radiance;
	}

	return result;
}

Vector LightSphere::calculate_intensity(Vector& light_vector, float pw)
{
	return radiance / pw;
}

Vector LightSphere::get_position(const Ray& ray, Intersection& intersection)
{
	Vector intersection_point = ray.origin + ray.direction * intersection.t;

	Vector local_intersection = inverse_transformation.transform_point(intersection_point);

	Vector d = rt_scene->vertex_data[center_vertex_id] - local_intersection;

	float distance = d.length();

	float sintheta_max = std::max(-1.0f, std::min(1.0f, radius / distance));

	float theta_max = asin(radius / distance);

	float costheta_max = cos(theta_max);

	float r1 = rt_scene->dist(rt_scene->rng);
	float r2 = rt_scene->dist(rt_scene->rng);

	costheta_max = std::max(0.0f, costheta_max);

	float phi = 2 * PI * r1;
	float theta = acos(1 - r2 + r2 * costheta_max);

	Vector w = d.normalized();

	Vector v = w.find_orthonormal_v();

	Vector u = w.cross(v).normalized();

	Vector l = w * cos(theta) + v * sin(theta) * cos(phi) + u * sin(theta) * sin(phi);

	l = base_transformation.transform_vector(l);

	Vector light_ray_origin = intersection_point + intersection.surface_normal * rt_scene->shadow_ray_epsilon;

	Ray light_ray = Ray(light_ray_origin, l);
	Intersection light_intersection;

	if (Sphere::intersect(light_ray, light_intersection, false))
	{
		intersection.pw = 1 / (2 * PI * (1 - costheta_max));
		return intersection_point + l * (light_intersection.t);
	}

	intersection.pw = 1e30;

	return intersection_point + l * light_intersection.t;
}