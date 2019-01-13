#pragma once
#ifndef __CAMERA__H_
#define __CAMERA__H_

#include "rtVector.h"
#include <string>
#include "ray.h"
#include "intersection.h"

class Camera
{
public:
	enum class Handedness
	{
		Right,
		Left
	};

	Camera();

	Vector sample_hemisphere(const Vector& direction, bool cosine_sampling);
	
	void initialize();
	void render();
	Vector calculate_pixel_color(Ray& ray);
	Vector send_ray(const Ray& ray, int level = 0, bool culling = true);
	bool send_shadow_ray(const Ray& ray, Intersection& intersection);
	Vector reflect_ray(const Ray& ray, const Intersection &intersection, int level, bool culling = true, bool inside = false);
	Vector refract_ray(const Ray& ray, const Intersection &intersection, int level);
	Vector calculate_lights(const Ray& ray, Intersection& intersection, int level);
	Vector trace_path(const Ray& ray, Intersection& intersection, int level);
	Vector calculate_diffuse_and_specular(const Ray& ray, Intersection& intersection);



	Vector position;
	Vector gaze;
	Vector up;
	Vector right;
	Vector near_plane;
	Vector top_left;

	float pixel_width;
	float pixel_height;
	float near_distance;
	float focus_distance;
	float aperture_size;
	int image_width;
	int image_height;

	int num_samples;

	bool is_hdr;
	std::string tmo;
	float image_key;
	float burn_percentage;
	float saturation;
	float gamma;

	std::string gamma_correction_method;

	std::string image_name;

	Handedness handedness;
};

#endif