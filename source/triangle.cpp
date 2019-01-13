#include "triangle.h"
#include "ray.h"
#include "scene.h"
#include <iostream>

Triangle::Triangle()
{
	v0_id = 0;
	v1_id = 0;
	v2_id = 0;
	vertex_offset = 0;
	texture_offset = 0;

	bbox.max = Vector();
	bbox.min = Vector();
}

Triangle::Triangle(int v0, int v1, int v2, int vertex_offset)
{
	dpdu = Vector(0.0f);
	dpdv = Vector(0.0f);

	v0_id = v0;
	v1_id = v1;
	v2_id = v2;

	Vector a = rt_scene->vertex_data[v0_id + vertex_offset];
	Vector b = rt_scene->vertex_data[v1_id + vertex_offset];
	Vector c = rt_scene->vertex_data[v2_id + vertex_offset];

	float minx = std::fmin(a.x, std::fmin(b.x, c.x));
	float miny = std::fmin(a.y, std::fmin(b.y, c.y));
	float minz = std::fmin(a.z, std::fmin(b.z, c.z));

	float maxx = std::fmax(a.x, std::fmax(b.x, c.x));
	float maxy = std::fmax(a.y, std::fmax(b.y, c.y));
	float maxz = std::fmax(a.z, std::fmax(b.z, c.z));

	bbox.min = Vector(minx, miny, minz);
	bbox.max = Vector(maxx, maxy, maxz);
}

void Triangle::calculate_dp()
{
	Vector p1 = rt_scene->vertex_data[v0_id + vertex_offset];
	Vector p2 = rt_scene->vertex_data[v1_id + vertex_offset];
	Vector p3 = rt_scene->vertex_data[v2_id + vertex_offset];

	Vector uv1 = rt_scene->texture_coordinate_data[v0_id + texture_offset];
	Vector uv2 = rt_scene->texture_coordinate_data[v1_id + texture_offset];
	Vector uv3 = rt_scene->texture_coordinate_data[v2_id + texture_offset];

	float a = uv2.x - uv1.x;
	float b = uv2.y - uv1.y;
	float c = uv3.x - uv1.x;
	float d = uv3.y - uv1.y;

	float coef = 1.0f / (a * d - b * c);

	a *= coef;
	b *= coef;
	c *= coef;
	d *= coef;

	Vector edge1 = p2 - p1;
	Vector edge2 = p3 - p1;

	dpdu = d * edge1 - b * edge2;
	dpdv = -c * edge1 + a * edge2;
}

bool Triangle::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	Vector direction = ray.direction;

	Vector a = rt_scene->vertex_data[v0_id + vertex_offset];
	Vector b = rt_scene->vertex_data[v1_id + vertex_offset];
	Vector c = rt_scene->vertex_data[v2_id + vertex_offset];

	Vector normal = ((b - a).cross(c - a)).normalized();

	if (culling && direction.dot(normal) > 0)
	{
		return false;
	}

	Vector edge1 = a - b;
	Vector edge2 = a - c;
	Vector replace_vector = a - ray.origin;

	float A = Vector::determinant_3(edge1, edge2, direction);

	float beta = Vector::determinant_3(replace_vector, edge2, direction) / A;
	if (beta < 0)
	{
		return false;
	}

	float gamma = Vector::determinant_3(edge1, replace_vector, direction) / A;
	if (gamma < 0 || beta + gamma > 1)
	{
		return false;
	}

	float t = Vector::determinant_3(edge1, edge2, replace_vector) / A;

	float alpha = 1.0f - beta - gamma;


	if (t > 1e-3 && intersection.t > t)
	{
		intersection.t = t;

		if (shading_mode == ShadingMode::Smooth)
		{
			normal = alpha * rt_scene->normal_data[v0_id + vertex_offset] + beta * rt_scene->normal_data[v1_id + vertex_offset] + gamma * rt_scene->normal_data[v2_id + vertex_offset];
			normal = normal.normalized();
		}

		intersection.surface_normal = normal;

		if (texture)
		{
			Vector u_vec = Vector(rt_scene->texture_coordinate_data[v0_id + texture_offset].x, rt_scene->texture_coordinate_data[v1_id + texture_offset].x, rt_scene->texture_coordinate_data[v2_id + texture_offset].x);
			Vector v_vec = Vector(rt_scene->texture_coordinate_data[v0_id + texture_offset].y, rt_scene->texture_coordinate_data[v1_id + texture_offset].y, rt_scene->texture_coordinate_data[v2_id + texture_offset].y);

			intersection.u = alpha * u_vec.x + beta * u_vec.y + gamma * u_vec.z;
			intersection.v = alpha * v_vec.x + beta * v_vec.y + gamma * v_vec.z;

			intersection.dpdu = dpdu;
			intersection.dpdv = dpdv;
		}

		return true;
	}

	return false;
}