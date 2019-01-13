#include "mesh.h"
#include <iostream>
#include <algorithm>
#include "scene.h"

Mesh::Mesh()
{
	bvh = nullptr;
	shading_mode = ShadingMode::Flat;
	texture_id = -1;

	motion_vector = Vector(0.0f);
	base_transformation = Matrix::Identity();
}

void Mesh::construct_bvh()
{
	BVH::SplitRange range = { 0, faces.size() };

	bvh = new BVH(faces, range);
	bbox = bvh->bbox;
}

bool Mesh::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	bool result = false;

	if (bvh->intersect(ray, intersection, culling, texture))
	{
		result = true;
		intersection.radiance = Vector(0.0f);
	}

	/*for (auto it = faces.begin(); it != faces.end(); it++)
	{
		if ((*it)->intersect(ray, intersection, culling))
		{
			result = true;
			intersection.radiance = Vector(0.0f);
		}
	}*/

	return result;
}

MeshInstance::MeshInstance()
{
	shading_mode = ShadingMode::Flat;
	texture_id = -1;

	motion_vector = Vector(0.0f);
	base_transformation = Matrix::Identity();
}

bool MeshInstance::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	bool result = false;

	Matrix inverse = inverse_transformation;
	Matrix transpose = transpose_transformation;

	if (motion_vector.x != 0.0f || motion_vector.y != 0.0f || motion_vector.z != 0.0f)
	{
		Matrix motion_matrix = Matrix::TranslationMatrix(motion_vector * ray.delta_time);

		Matrix transformation = motion_matrix * base_transformation;
		inverse = transformation.inverse();
		transpose = inverse.transpose();
	}

	Ray transformed_ray = Ray(inverse.transform_point(ray.origin), inverse.transform_vector(ray.direction));
	transformed_ray.inside = ray.inside;

	bool need_uv = (texture_id != -1) && (rt_scene->textures[texture_id].texture_type == Texture::TextureType::Image);

	if (rt_scene->meshes[base_mesh_id]->intersect(transformed_ray, intersection, culling, need_uv))
	{
		result = true;
		intersection.material_id = material_id;

		if(texture_id != -1)
		{
			Vector intersection_point = transformed_ray.origin + transformed_ray.direction * intersection.t;
			if (rt_scene->textures[texture_id].bump_map)
			{
				Vector dpdu = intersection.dpdu;
				Vector dpdv = intersection.dpdv;

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
	}

	return result;
}

//Light Mesh
LightMesh::LightMesh()
{
	bvh = nullptr;
	shading_mode = ShadingMode::Flat;
	texture_id = -1;

	motion_vector = Vector(0.0f);
	base_transformation = Matrix::Identity();

}

void LightMesh::calculate_area()
{
	total_area = 0;

	for (int i = 0; i < faces.size(); i++)
	{
		Triangle *face = (Triangle*)faces[i];

		Vector p1 = rt_scene->vertex_data[face->v0_id + face->vertex_offset];
		Vector p2 = rt_scene->vertex_data[face->v1_id + face->vertex_offset];
		Vector p3 = rt_scene->vertex_data[face->v2_id + face->vertex_offset];

		p1 = base_transformation.transform_point(p1);
		p2 = base_transformation.transform_point(p2);
		p3 = base_transformation.transform_point(p3);

		Vector edge1 = p2 - p1;
		Vector edge2 = p3 - p1;

		float area = (edge1.cross(edge2)).length() / 2;

		total_area += area;

		areas[i] += total_area;
	}
}

bool LightMesh::intersect(const Ray& ray, Intersection& intersection, bool culling, bool texture)
{
	bool result = false;

	if (bvh->intersect(ray, intersection, culling, texture))
	{
		result = true;
		intersection.radiance = radiance;
	}

	return result;
}

Vector LightMesh::calculate_intensity(Vector& light_vector, float pw)
{

	//std::cout << 2 << pw << std::endl;

	return radiance / pw;
}

Vector LightMesh::get_position(const Ray& ray, Intersection& intersection)
{
	Triangle *face = (Triangle*)faces[0];

	float cdf = rt_scene->dist(rt_scene->rng) * total_area;

	for (int i = 0; i < faces.size(); i++)
	{
		if (cdf < areas[i])
		{
			face = (Triangle*)faces[i];
			break;
		}
	}

	Vector p1 = rt_scene->vertex_data[face->v0_id + face->vertex_offset];
	Vector p2 = rt_scene->vertex_data[face->v1_id + face->vertex_offset];
	Vector p3 = rt_scene->vertex_data[face->v2_id + face->vertex_offset];

	p1 = base_transformation.transform_point(p1);
	p2 = base_transformation.transform_point(p2);
	p3 = base_transformation.transform_point(p3);

	float r1 = rt_scene->dist(rt_scene->rng);
	float r2 = sqrt(rt_scene->dist(rt_scene->rng));

	Vector q = (1 - r1) * p2 + r1 * p3;
	Vector p = (1 - r2) * p1 + r2 * q;

	Vector edge1 = p2 - p1;
	Vector edge2 = p3 - p1;

	Vector normal = edge1.cross(edge2).normalized();

	Vector position = p + normal * rt_scene->shadow_ray_epsilon;

	position = base_transformation.transform_point(position);

	Vector intersection_point = (ray.origin + ray.direction * intersection.t);

	Vector light = intersection_point - position;

	float distance_squared = light.length_squared();

	light = light.normalized();

	float cosine = std::max(0.0f, normal.dot(light));

	intersection.pw = distance_squared / (cosine * total_area);

	return position;
}