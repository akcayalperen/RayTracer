#pragma once
#ifndef __TEXTURE__H_
#define __TEXTURE__H_

#include "rtVector.h"
#include <string>

class Texture
{
public:
	enum class DecalMode
	{
		None,
		ReplaceKD,
		BlendKD,
		ReplaceAll
	};

	enum class Appearance
	{
		Clamp,
		Repeat,
		Patch,
		Vein
	};

	enum class TextureType
	{
		Image,
		Perlin
	};

	Texture();

	Vector get_color(float u, float v);
	Vector get_perlin_color(float x, float y, float z);
	Vector get_bump(float u, float v);
	Vector get_perlin_bump(float x, float y, float z);
	Vector bilinear_interpolation(float u, float v, bool repeat);
	Vector nearest_neighbor(float u, float v);
	void set_decal_mode(std::string mode);
	void set_appearance(std::string mode);


	std::string image_name;
	std::string interpolation;
	std::string decal_mode_string;
	std::string appearance_string;
	DecalMode decal_mode;
	Appearance apperance;
	float normalizer;

	float scaling_factor;
	TextureType texture_type;

	unsigned char *image_data;
	float *radiance_data;
	int width;
	int height;
	int number_of_components;

	bool bump_map;
	float bump_map_multiplier;

	bool degamma;
};

#endif