#pragma once
#ifndef __SCENE__H_
#define __SCENE__H_

#include <vector>
#include "rtMatrix.h"
#include "rtVector.h"
#include "camera.h"
#include "light.h"
#include "material.h"
#include "texture.h"
#include "object.h"
#include "BRDF.h"
#include <random>

#define STB_IMAGE_IMPLEMENTATION

class Scene
{
public:
	enum class Integrator
	{
		RayTracing,
		PathTracing
	};

	Scene();
	void construct_bvh();
	void render();
	Vector get_background_color(float i, float j);

	Object *bvh;

	std::mt19937 rng;
	std::uniform_real_distribution<double> dist;

	Texture* background_texture;
	Texture* environment_map;

	Vector background_color;
	int max_recursion_depth;
	float shadow_ray_epsilon;
	float intersection_test_epsilon;

	std::vector<Camera> cameras;
	Vector ambient_light;

	std::vector<Light*> lights;
	std::vector<Material> materials;
	std::vector<Texture> textures;
	std::vector<Vector> vertex_data;
	std::vector<Vector> texture_coordinate_data;
	std::vector<Vector> normal_data;

	std::vector<BRDF*> brdfs;

	std::vector<Object*> objects;
	std::vector<Object*> meshes;
	std::vector<Object*> spheres;

	std::vector<Matrix> translations;
	std::vector<Matrix> scalings;
	std::vector<Matrix> rotations;

	Integrator integrator;
	bool importance_sampling;

	bool zero_index;
};

extern Scene* rt_scene;

#endif //__SCENE__H_

