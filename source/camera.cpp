#include "camera.h"
#include "scene.h"
#include "lodepng.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "ctpl_stl.h"
#include <time.h>
#include <stdlib.h>

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

#define THREAD_COUNT 16

inline float gaussian(float x, float y, float sigma)
{
	float r = std::sqrt(x*x + y*y);
	float s = 2.0f * sigma * sigma;
	return (std::exp(-(r*r) / s)) / PI * s;
}

inline float clamp(const float &lo, const float &hi, const float &v)
{
	return std::max(lo, std::min(hi, v));
}

Camera::Camera()
{
	position = Vector(true);
	gaze = Vector();
	up = Vector();
	right = Vector();
	near_plane = Vector();

	focus_distance = -1.0f;
	aperture_size = 0;

	num_samples = 9;

	is_hdr = false;
	gamma = 2.2f;

	gamma_correction_method = "";

	handedness = Handedness::Right;
}

void Camera::initialize()
{
	gaze = gaze.normalized();
	up = up.normalized();

	
	if (handedness == Handedness::Right)
	{
		right = gaze.cross(up).normalized();
		up = right.cross(gaze).normalized();
	}
	else
	{
		right = up.cross(gaze).normalized();
		up = gaze.cross(right).normalized();
	}

	if (focus_distance > 0.0f)
	{
		float focal_ratio = focus_distance / near_distance;
		near_distance = focus_distance;
		
		near_plane.x *= focal_ratio;
		near_plane.y *= focal_ratio;
		near_plane.z *= focal_ratio;
		near_plane.w *= focal_ratio;
	}

	top_left = position + near_distance * gaze + near_plane.x * right + near_plane.w * up;
	pixel_width = (near_plane.y - near_plane.x) / (float)image_width;
	pixel_height = (near_plane.w - near_plane.z) / (float)image_height;
}

Vector UniformSampleHemisphere(float u1, float u2)
{
	const float r = sqrtf(1.0f - u1 * u1);
	const float phi = 2 * PI * u2;

	return Vector(cos(phi) * r, sin(phi) * r, u1);
}

Vector Camera::sample_hemisphere(const Vector& direction, bool cosine_sampling)
{
	float r1 = rt_scene->dist(rt_scene->rng);
	float r2 = rt_scene->dist(rt_scene->rng);

	float phi;
	float theta;

	if (cosine_sampling)
	{
		phi = 2 * PI * r1;
		theta = asin(sqrt(r2));
	}
	else
	{
		phi = 2 * PI * r1;
		theta = acos(r2);
	}

	Vector w = direction;
	w = w.normalized();

	Vector v = w.find_orthonormal_v();

	Vector u = w.cross(v).normalized();

	Vector l = w * cos(theta) + v * sin(theta) * cos(phi) + u * sin(theta) * sin(phi);

	l = l.normalized();

	return l;
}

Vector Camera::calculate_pixel_color(Ray& ray)
{
	ray.delta_time = rt_scene->dist(rt_scene->rng);// ((double)rand() / RAND_MAX);

	return send_ray(ray);
}

Vector Camera::send_ray(const Ray& ray, int level, bool culling)
{
	bool result = false;
	Intersection intersection;

	if (rt_scene->bvh && rt_scene->bvh->intersect(ray, intersection, culling))
	{
		result = true;
	}

	for (auto it = rt_scene->spheres.begin(); it != rt_scene->spheres.end(); it++)
	{
		if ((*it)->intersect(ray, intersection, culling))
		{
			result = true;
		}
	}

	if (result)
	{
		if (level >= rt_scene->max_recursion_depth)
		{
			return intersection.radiance;
		}

		if (rt_scene->integrator == Scene::Integrator::PathTracing)
		{
			return trace_path(ray, intersection, level);
		}
		else
		{
			return calculate_lights(ray, intersection, level);
		}
	}

	if (level == 0 && rt_scene->integrator == Scene::Integrator::RayTracing)
	{
		if (rt_scene->environment_map)
		{
			Vector l = ray.direction;
			l = l.normalized();

			float theta = acos(l.y);
			float phi = atan2(l.z, l.x);

			float u_coord = (-phi + PI) / (2 * PI);
			float v_coord = theta / PI;

			return rt_scene->environment_map->get_color(u_coord, v_coord);

		}
		else
		{
			return rt_scene->get_background_color(ray.target_pixel_x, ray.target_pixel_y);
		}
	}

	return Vector(0.0f);
}

bool Camera::send_shadow_ray(const Ray& ray, Intersection& intersection)
{
	bool result = false;

	if (rt_scene->bvh && rt_scene->bvh->intersect(ray, intersection, false))
	{
		result = true;
	}

	for (auto it = rt_scene->spheres.begin(); it != rt_scene->spheres.end(); it++)
	{
		if ((*it)->intersect(ray, intersection, false))
		{
			result = true;
		}
	}

	return result;
}

Vector Camera::reflect_ray(const Ray& ray, const Intersection& intersection, int level, bool culling, bool inside)
{
	Vector normalized_direction = ray.direction;
	normalized_direction = normalized_direction.normalized();

	Vector glossiness = Vector(0.0f);

	float roughness = rt_scene->materials[intersection.material_id].roughness;

	Vector surface_normal = intersection.surface_normal;
	Vector w_r = (normalized_direction + intersection.surface_normal * 2 * surface_normal.dot(-normalized_direction)).normalized();

	if (roughness > 0.0f)
	{
		Vector rI;

		if (w_r.x < w_r.y && w_r.x < w_r.y)
		{
			rI = Vector(1.0f, w_r.y, w_r.z);
		}
		else if (w_r.y < w_r.x && w_r.y < w_r.z)
		{
			rI = Vector(w_r.x, 1.0f, w_r.z);
		}
		else
		{
			rI = Vector(w_r.x, w_r.y, 1.0f);
		}

		double r1 = rt_scene->dist(rt_scene->rng);//((double)rand() / RAND_MAX);
		double r2 = rt_scene->dist(rt_scene->rng);// ((double)rand() / RAND_MAX);

		Vector u = rI.cross(w_r);
		Vector v = u.cross(w_r);

		glossiness = (u*r1 + v*r2)*roughness;
	}
	
	w_r = (w_r + glossiness).normalized();

	Ray mirror_ray = Ray(ray.origin + ray.direction * intersection.t + w_r * rt_scene->shadow_ray_epsilon, w_r);
	mirror_ray.inside = inside;
	mirror_ray.delta_time = ray.delta_time;
	mirror_ray.light_hit = ray.light_hit;

	return send_ray(mirror_ray, level + 1, culling);
}

Vector Camera::refract_ray(const Ray& ray, const Intersection& intersection, int level)
{
	Vector transparency = rt_scene->materials[intersection.material_id].transparency;
	float refraction_index = rt_scene->materials[intersection.material_id].refraction_index;

	Vector intersection_point = ray.origin + ray.direction * intersection.t;

	Vector direction = ray.direction;
	direction = direction.normalized();

	Vector reflection_color = Vector(0.0f);
	Vector refraction_color = Vector(0.0f);

	Vector normal = intersection.surface_normal;

	float costheta;
	float c;

	Vector k = Vector(1.0f);

	if (direction.dot(normal) < 0)
	{
		reflection_color = reflect_ray(ray, intersection, level, true, false);

		costheta = (-direction).dot(normal);

		float index_ratio = 1.0f / refraction_index;
		float delta = 1.0f - index_ratio * index_ratio * (1.0f - costheta * costheta);

		Vector refraction_direction = (direction + normal * costheta) * index_ratio - normal * sqrtf(delta);
		refraction_direction = refraction_direction.normalized();
		Ray refraction_ray = Ray(intersection_point + refraction_direction * rt_scene->shadow_ray_epsilon, refraction_direction);
		//refraction_ray.inside = true;
		refraction_ray.delta_time = ray.delta_time;

		refraction_color = send_ray(refraction_ray, level + 1, false);
		c = costheta;
	}
	else
	{
		reflection_color = reflect_ray(ray, intersection, level, false, true);

		normal = -normal;
		costheta = (-direction).dot(normal);
		k.x = std::exp(std::log(transparency.x) * intersection.t);
		k.y = std::exp(std::log(transparency.y) * intersection.t);
		k.z = std::exp(std::log(transparency.z) * intersection.t);

		float index_ratio = refraction_index;
		float delta = 1.0f - index_ratio * index_ratio * (1.0f - costheta * costheta);
		if (delta < 0.0f)
		{
			return k * reflection_color;
		}

		Vector refraction_direction = (direction + normal * costheta) * index_ratio - normal * sqrtf(delta);
		refraction_direction = refraction_direction.normalized();
		Ray refraction_ray = Ray(intersection_point + refraction_direction * rt_scene->shadow_ray_epsilon, refraction_direction);
		//refraction_ray.inside = false;
		refraction_ray.delta_time = ray.delta_time;
		refraction_ray.light_hit = ray.light_hit;
		
		refraction_color = send_ray(refraction_ray, level + 1, true);
		c = refraction_direction.dot(-normal);
	}

	float R0 = (refraction_index - 1.0f) / (refraction_index + 1.0f);
	R0 = R0 * R0;
	float R = R0 + (1.0f - R0) * std::pow((1.0f - c), 5);

	return k*(R * reflection_color + (1.0f - R) * refraction_color);
}

Vector Camera::calculate_lights(const Ray& ray, Intersection& intersection, int level)
{
	Vector mirror_intensity = Vector(0.0f);
	Vector refraction_intensity = Vector(0.0f);
	Vector ambient_intensity = Vector(0.0f); 
	Vector diffuse_specular_intensity = Vector(0.0f);

	if (!ray.inside)
	{
		ambient_intensity = rt_scene->ambient_light * rt_scene->materials[intersection.material_id].ambient;
	
		diffuse_specular_intensity = calculate_diffuse_and_specular(ray, intersection);
	}

	Vector mirror_effect = rt_scene->materials[intersection.material_id].mirror;
	Vector transparency = rt_scene->materials[intersection.material_id].transparency;
	if (level < rt_scene->max_recursion_depth)
	{
		if (mirror_effect.x + mirror_effect.y + mirror_effect.z > 0.0f)
		{
			mirror_intensity = mirror_effect * reflect_ray(ray, intersection, level);
		}
		if (transparency.x + transparency.y + transparency.z > 0.0f)
		{
			refraction_intensity = refract_ray(ray, intersection, level);
		}
	}

	return ambient_intensity + diffuse_specular_intensity + mirror_intensity + refraction_intensity;
}

Vector Camera::trace_path(const Ray& ray, Intersection& intersection, int level)
{
	Vector radiance = Vector(0.0f);

	if (intersection.radiance.x > 0.0f)
	{
		if (ray.light_hit)
		{
			return Vector(0.0f);
		}
		return intersection.radiance;
	}

	Vector direct_contribution = calculate_diffuse_and_specular(ray, intersection) / (rt_scene->lights.size());

	Ray new_ray = ray;

	if (direct_contribution.x > 0.0f || direct_contribution.y > 0.0f || direct_contribution.z > 0.0f)
	{
		new_ray.light_hit = true;
	}

	radiance = radiance + direct_contribution;

	Material material = rt_scene->materials[intersection.material_id];

	if (!ray.inside)
	{
		Vector intersection_point = ray.origin + ray.direction * intersection.t;
		Vector l = sample_hemisphere(intersection.surface_normal, rt_scene->importance_sampling);
		Ray diffuse_ray = Ray(intersection_point + intersection.surface_normal * rt_scene->shadow_ray_epsilon * 0.1f, l);
		diffuse_ray.inside = ray.inside;
		diffuse_ray.light_hit = ray.light_hit;

		float cosinetheta = std::max(0.0f, l.dot(intersection.surface_normal));

		Vector normalized_direction = ray.direction;
		normalized_direction = normalized_direction.normalized();

		Vector half_vector = (l - normalized_direction).normalized();
		float cosinealpha = std::max(0.0f, half_vector.dot(intersection.surface_normal));

		Vector diffuse = material.diffuse;
		Vector specular = material.specular;

		diffuse = diffuse / PI;

		float exponent = 50;
		specular = specular * ((exponent + 8) / (8 * PI));

		if (rt_scene->importance_sampling)
		{
			if (intersection.decal_mode == Texture::DecalMode::ReplaceKD)
			{
				diffuse = intersection.texture_color;
				specular = material.specular;
				exponent = material.phong_exponent;
			}

			Vector diffuse_specular = (diffuse + specular *
				pow(cosinealpha, exponent)) * send_ray(diffuse_ray, level + 1) * PI;

			radiance = radiance + diffuse_specular;
		}
		else
		{
			Vector diffuse_specular = (diffuse + specular *
				pow(cosinealpha, exponent)) * cosinetheta * send_ray(diffuse_ray, level + 1) * 2 * PI;

			radiance = radiance + diffuse_specular;
		}
	}

	Vector mirror_effect = material.mirror;
	Vector transparency = material.transparency;
	if (level < rt_scene->max_recursion_depth)
	{
		if (mirror_effect.x + mirror_effect.y + mirror_effect.z > 0.0f)
		{
			radiance = radiance + mirror_effect * reflect_ray(new_ray, intersection, level);
		}
		if (transparency.x + transparency.y + transparency.z > 0.0f)
		{
			radiance = radiance + refract_ray(new_ray, intersection, level);
		}
	}

	

	return radiance;
}

Vector Camera::calculate_diffuse_and_specular(const Ray& ray, Intersection& intersection)
{
	Vector diffuse_color = Vector(0.0f);
	Vector specular_color = Vector(0.0f);
	Vector intersection_point = ray.origin + ray.direction * intersection.t;

	if (intersection.radiance.x > 0.0f || intersection.radiance.y > 0.0f || intersection.radiance.z > 0.0f)
	{
		return intersection.radiance;
	}

	if (rt_scene->environment_map)
	{
		Vector l = sample_hemisphere(intersection.surface_normal, rt_scene->importance_sampling);

		float theta = acos(l.y);
		float phi = atan2(l.z, l.x);

		float u_coord = (-phi + PI) / (2 * PI);
		float v_coord = theta / PI;

		float cosine = std::max(0.0f, l.dot(intersection.surface_normal));


		Vector shadow_ray_origin = intersection_point + (intersection.surface_normal * rt_scene->shadow_ray_epsilon);
		Ray shadow_ray = Ray(shadow_ray_origin, l);
		shadow_ray.delta_time = ray.delta_time;
		Intersection shadow_intersection;

		if (send_shadow_ray(shadow_ray, shadow_intersection))
		{
			diffuse_color = Vector(0.0f);
		}
		else
		{
			Vector color = intersection.texture_color * cosine * rt_scene->environment_map->get_color(u_coord, v_coord);

			diffuse_color = diffuse_color + color * (2 * PI);
		}
	}

	for (auto light = rt_scene->lights.begin(); light != rt_scene->lights.end(); light++)
	{
		if (intersection.decal_mode == Texture::DecalMode::ReplaceAll)
		{
			return intersection.texture_color * 255.0f;
		}

		Vector light_position = (*light)->get_position(ray, intersection);

		Vector light_vector = light_position - intersection_point;
		Vector light_direction = light_vector.normalized();

		auto direction = ray.direction;

		Vector half_vector = (light_direction - direction.normalized()).normalized();
		float light_distance = light_vector.length();

		//Is Facing Light
		Vector surface_normal = intersection.surface_normal;
		if (surface_normal.dot(light_direction) >= 0.0f)
		{
			Vector shadow_ray_origin = intersection_point + (intersection.surface_normal * rt_scene->shadow_ray_epsilon);
			Ray shadow_ray = Ray(shadow_ray_origin, light_direction);
			shadow_ray.delta_time = ray.delta_time;
			Intersection shadow_intersection;

			if (send_shadow_ray(shadow_ray, shadow_intersection))
			{
				float shadow_distance = (shadow_ray.direction * shadow_intersection.t).length();
				if (shadow_distance < light_distance - rt_scene->shadow_ray_epsilon)
				{
					continue;
				}
			}

			int BRDF_id = rt_scene->materials[intersection.material_id].brdf_id;

			float cosine = std::max(0.0f, light_direction.dot(surface_normal));

			if (BRDF_id != -1)
			{
				diffuse_color = diffuse_color + rt_scene->brdfs[BRDF_id]->shade(ray, intersection, *light, light_vector);
			}
			else
			{
				if (intersection.decal_mode == Texture::DecalMode::ReplaceKD)
				{
					diffuse_color = diffuse_color + intersection.texture_color * cosine * (*light)->calculate_intensity(light_vector, cosine);
				}
				else if (intersection.decal_mode == Texture::DecalMode::BlendKD)
				{
					Vector blend = (intersection.texture_color + rt_scene->materials[intersection.material_id].diffuse) / 2;
					diffuse_color = diffuse_color + blend * cosine * (*light)->calculate_intensity(light_vector, cosine);
				}
				else
				{
					diffuse_color = diffuse_color + rt_scene->materials[intersection.material_id].diffuse * cosine * (*light)->calculate_intensity(light_vector, cosine);
				}

				cosine = std::max(0.0f, half_vector.dot(surface_normal));
				specular_color = specular_color + rt_scene->materials[intersection.material_id].specular *
					pow(cosine, rt_scene->materials[intersection.material_id].phong_exponent) *
					(*light)->calculate_intensity(light_vector, cosine);
			}
		}
	}

	return diffuse_color + specular_color;
}

void Camera::render()
{
	const int width = image_width;
	const int height = image_height;

	const int samples = (int)std::sqrt(num_samples);
	const float sigma = 1.0f / samples;

	unsigned char *image = new unsigned char[width * height * 3];
	float *hdr_image = new float[width * height * 3];

	//int y = 400;
	//int x = 400;
	//srand(time(NULL));

	ctpl::thread_pool pool(THREAD_COUNT);
	for (int y = 0; y < height; y++)
	{
		pool.push([y, &width, &height, &samples, &sigma, &image, &hdr_image, this](int id)
		{
			for (int x = 0; x < width; x++)
			{

				Vector color = Vector();
				float total_weight = 0.0f;

				for (int i = 0; i < samples; i++)
				{
					for (int j = 0; j < samples; j++)
					{
						float r1 = (i + /*((double)rand() / RAND_MAX)*/ rt_scene->dist(rt_scene->rng)) / (float)samples - 0.5f;
						float r2 = (j + /*((double)rand() / RAND_MAX)*/ rt_scene->dist(rt_scene->rng)) / (float)samples - 0.5f;

						Vector horizontal_aperture = Vector(0.0f);
						Vector vertical_aperture = Vector(0.0f);
						if (aperture_size > 0.0f)
						{
							float r3 = (/*(double)rand() / RAND_MAX*/rt_scene->dist(rt_scene->rng)) - (aperture_size / 2);
							float r4 = (/*(double)rand() / RAND_MAX*/rt_scene->dist(rt_scene->rng)) - (aperture_size / 2);

							horizontal_aperture = r3 * right;
							vertical_aperture = r4 * up;
						}

						Vector destination = this->top_left + this->pixel_width * (0.5f + x + r1) * right - pixel_height * (0.5f + y + r2) * up;

						Vector ray_origin = position + vertical_aperture + horizontal_aperture;
						Ray primaryRay = Ray(ray_origin, destination - ray_origin);
						primaryRay.target_pixel_x = (float)x / width;
						primaryRay.target_pixel_y = (float)y / height;

						//float weight = gaussian(r1, r2, sigma);
						color = color + calculate_pixel_color(primaryRay);// *weight;
						//total_weight += weight;
					}
				}
				color = color / (float)num_samples;//total_weight;

				int index = (y * width + x) * 3;

				if (is_hdr)
				{
					hdr_image[index] = color.x;
					hdr_image[index + 1] = color.y;
					hdr_image[index + 2] = color.z;
				}

				color.clamp(0, 255.f);
				image[index] = (unsigned char)color.x;
				image[index + 1] = (unsigned char)color.y;
				image[index + 2] = (unsigned char)color.z;
			}
		});
	}
	pool.stop(true);

	if (is_hdr)
	{
		//Move to Function
		std::vector<float> luminances;

		float average_luminance = 0;
		int total_pixels = width*height;

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int index = (y * width + x) * 3;

				float luminance = 0.21 * hdr_image[index] + 0.72 * hdr_image[index + 1] + 0.07 * hdr_image[index + 2];

				luminances.push_back(luminance);

				average_luminance += std::log(luminance + 0.001);
			}
		}

		float log_average_luminance = std::exp(average_luminance / total_pixels);

		for (int i = 0; i < luminances.size(); i++)
		{
			luminances[i] = image_key * luminances[i] / log_average_luminance;
		}

		std::sort(luminances.begin(), luminances.end());

		int l_white_index = (int)((burn_percentage / 100.0f) * total_pixels);
		float l_white = luminances[total_pixels - l_white_index - 1];
		float l_white_squared = l_white * l_white;

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				int index = (y * width + x) * 3;

				float luminance = 0.21 * hdr_image[index] + 0.72 * hdr_image[index + 1] + 0.07 * hdr_image[index + 2];

				float luminance_scaled = luminance * (image_key / (log_average_luminance));

				float l_final = (luminance_scaled * (1 + (luminance_scaled / l_white_squared))) / (1 + luminance_scaled);

				float R = 0;
				float G = 0;
				float B = 0;
				if (luminance != 0)
				{
					R = std::pow((float)hdr_image[index] / luminance, saturation) * l_final;
					G = std::pow((float)hdr_image[index + 1] / luminance, saturation) * l_final;
					B = std::pow((float)hdr_image[index + 2] / luminance, saturation) * l_final;
				}

				if (gamma_correction_method == "sRGB")
				{
					image[index] = clamp(0.0f, 1.0f, (1.055f * std::pow(clamp(0.0f, 1.0f, R), 1 / 2.4f) - 0.055f)) * 255;
					image[index + 1] = clamp(0.0f, 1.0f, (1.055f * std::pow(clamp(0.0f, 1.0f, G), 1 / 2.4f) - 0.055f)) * 255;
					image[index + 2] = clamp(0.0f, 1.0f, (1.055f * std::pow(clamp(0.0f, 1.0f, B), 1 / 2.4f) - 0.055f)) * 255;
				}
				else
				{
					image[index] = std::pow(clamp(0.0f, 1.0f, R), 1 / gamma) * 255;
					image[index + 1] = std::pow(clamp(0.0f, 1.0f, G), 1 / gamma) * 255;
					image[index + 2] = std::pow(clamp(0.0f, 1.0f, B), 1 / gamma) * 255;
				}
			}
		}

		//Write EXR
		{
			EXRHeader header;
			InitEXRHeader(&header);

			EXRImage exr_image;
			InitEXRImage(&exr_image);

			exr_image.num_channels = 3;

			std::vector<float> images[3];
			images[0].resize(width * height);
			images[1].resize(width * height);
			images[2].resize(width * height);

			for (int i = 0; i < width * height; i++)
			{
				images[0][i] = hdr_image[3 * i + 0];
				images[1][i] = hdr_image[3 * i + 1];
				images[2][i] = hdr_image[3 * i + 2];
			}

			float* image_ptr[3];
			image_ptr[0] = &(images[2].at(0)); // B
			image_ptr[1] = &(images[1].at(0)); // G
			image_ptr[2] = &(images[0].at(0)); // R

			exr_image.images = (unsigned char**)image_ptr;
			exr_image.width = width;
			exr_image.height = height;

			header.num_channels = 3;
			header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
			// Must be (A)BGR order, since most of EXR viewers expect this channel order.
			strncpy(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
			strncpy(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
			strncpy(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';

			header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
			header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
			for (int i = 0; i < header.num_channels; i++) {
				header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
				header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
			}

			std::string exr_image_name = image_name + ".exr";

			const char* err;
			int ret = SaveEXRImageToFile(&exr_image, &header, exr_image_name.c_str(), &err);
			if (ret != TINYEXR_SUCCESS) {
				fprintf(stderr, "Save EXR err: %s\n", err);
			}

			printf("Saved exr file. [ %s ] \n", exr_image_name.c_str());

			free(header.channels);
			free(header.pixel_types);
			free(header.requested_pixel_types);
		}
	}

	lodepng::encode(image_name + ".png", image, width, height, LCT_RGB);
}