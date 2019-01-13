#include "scene.h"
#include "BVH.h"

Scene::Scene()
{
	background_color = Vector();
	max_recursion_depth = 0;
	shadow_ray_epsilon = 0.001f;
	intersection_test_epsilon = 0.001f;

	rng.seed(std::random_device()());
	dist = std::uniform_real_distribution<double>(0.0, 1.0);

	ambient_light = Vector();
	bvh = nullptr;
	background_texture = nullptr;
	environment_map = nullptr;

	integrator == Integrator::RayTracing;
	importance_sampling = false;

	zero_index = false;
}

void Scene::construct_bvh()
{
	BVH::SplitRange range = { 0, objects.size() };

	bvh = new BVH(objects, range);
}

void Scene::render()
{
	if(objects.size() > 0)
	{
		construct_bvh();
	}

	for (auto it = cameras.begin(); it != cameras.end(); it++)
	{
		it->initialize();

		it->render();
	}
}

Vector Scene::get_background_color(float i, float j)
{
	if (background_texture != nullptr)
	{
		return background_texture->get_color(i, j);
	}
	
	return background_color;
}