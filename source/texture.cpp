#include "texture.h"
#include <algorithm>
#include <iostream>

//TODO: Move This Table to Elsewhere
Vector gradient_table[16] = {
	Vector(1.0f, 1.0f, 0.0f),  Vector(-1.0f, 1.0f ,0.0f),
	Vector(1.0f, -1.0f, 0.0f), Vector(-1.0f, -1.0f, 0.0f),
	Vector(1.0f, 0.0f, 1.0f),  Vector(-1.0f, 0.0f, 1.0f),
	Vector(1.0f, 0.0f, -1.0f), Vector(-1.0f, 0.0f, -1.0f),
	Vector(0.0f, 1.0f, 1.0f),  Vector(0.0f, -1.0f, 1.0f),
	Vector(0.0f, 1.0f, -1.0f), Vector(0.0f, -1.0f, -1.0f),
	Vector(1.0f, 1.0f, 0.0f),  Vector(-1.0f, 1.0f, 0.0f),
	Vector(0.0f, -1.0f, 1.0f), Vector(0.0f, -1.0f, -1.0f)
};

int index_table[16] = { 1, 6, 5, 13, 11, 0, 8, 2, 4, 10, 9, 14, 15, 12, 7, 3 };

static inline float fade(float x)
{
	x = abs(x);
	return -6 * pow(x, 5) + 15 * pow(x, 4) - 10 * pow(x, 3) + 1;
}

static inline int permutation(int i)
{
	return index_table[(16 + (i % 16)) % 16];
}

Texture::Texture()
{
	normalizer = 255.0f;
	image_data = nullptr;
	radiance_data = nullptr;
	decal_mode = DecalMode::None;
	apperance = Appearance::Repeat;
	scaling_factor = 1.0f;
	bump_map = false;
	bump_map_multiplier = 0.0f;
	degamma = false;
}

Vector Texture::nearest_neighbor(float i, float j)
{
	int x = round(i);
	int y = round(j);
	if (x == height)
	{
		x--;
	}
	if (y == width)
	{
		y--;
	}
	int index = ((y * width + x) * number_of_components);
	return Vector(image_data[index], image_data[index + 1], image_data[index + 2]);
}

Vector Texture::bilinear_interpolation(float i, float j, bool repeat)
{
	int p = floor(i);
	int q = floor(j);
	float dx = i - p;
	float dy = j - q;

	int top_left_pixel = ((q % height) * width + (p % width)) * number_of_components;
	int top_right_pixel = ((q % height) * width + ((p + 1) % width)) * number_of_components;
	int bottom_left_pixel = (((q + 1) % height) * width + (p % width)) * number_of_components;
	int bottom_right_pixel = (((q + 1) % height) * width + ((p + 1) % width)) * number_of_components;

	if (image_data != nullptr)
	{
		Vector top_left_color(image_data[top_left_pixel], image_data[top_left_pixel + 1], image_data[top_left_pixel + 2]);
		Vector top_right_color(image_data[top_right_pixel], image_data[top_right_pixel + 1], image_data[top_right_pixel + 2]);
		Vector bottom_left_color(image_data[bottom_left_pixel], image_data[bottom_left_pixel + 1], image_data[bottom_left_pixel + 2]);
		Vector bottom_right_color(image_data[bottom_right_pixel], image_data[bottom_right_pixel + 1], image_data[bottom_right_pixel + 2]);

		return top_left_color * (1 - dx) * (1 - dy) +
			top_right_color * (dx) * (1 - dy) +
			bottom_left_color * (1 - dx) * (dy)+
			bottom_right_color * (dx) * (dy);
	}
	else if (radiance_data != nullptr)
	{
		Vector top_left_color(radiance_data[top_left_pixel], radiance_data[top_left_pixel + 1], radiance_data[top_left_pixel + 2]);
		Vector top_right_color(radiance_data[top_right_pixel], radiance_data[top_right_pixel + 1], radiance_data[top_right_pixel + 2]);
		Vector bottom_left_color(radiance_data[bottom_left_pixel], radiance_data[bottom_left_pixel + 1], radiance_data[bottom_left_pixel + 2]);
		Vector bottom_right_color(radiance_data[bottom_right_pixel], radiance_data[bottom_right_pixel + 1], radiance_data[bottom_right_pixel + 2]);

		return top_left_color * (1 - dx) * (1 - dy) +
			top_right_color * (dx) * (1 - dy) +
			bottom_left_color * (1 - dx) * (dy)+
			bottom_right_color * (dx) * (dy);
	}
	return Vector(0.0f);
}

Vector Texture::get_perlin_color(float x, float y, float z)
{
	x *= scaling_factor;
	y *= scaling_factor;
	z *= scaling_factor;

	int i = (int)floor(x);
	int j = (int)floor(y);
	int k = (int)floor(z);

	Vector edge1 = gradient_table[permutation(i + permutation(j + permutation(k)))];
	Vector edge2 = gradient_table[permutation(i + 1 + permutation(j + permutation(k)))];
	Vector edge3 = gradient_table[permutation(i + 1 + permutation(j + 1 + permutation(k)))];
	Vector edge4 = gradient_table[permutation(i + permutation(j + 1 + permutation(k)))];
	Vector edge5 = gradient_table[permutation(i + permutation(j + permutation(k + 1)))];
	Vector edge6 = gradient_table[permutation(i + 1 + permutation(j + permutation(k + 1)))];
	Vector edge7 = gradient_table[permutation(i + 1 + permutation(j + 1 + permutation(k + 1)))];
	Vector edge8 = gradient_table[permutation(i + permutation(j + 1 + permutation(k + 1)))];

	Vector v1 = Vector(x - i, y - j, z - k);
	Vector v2 = Vector(x - (i + 1), y - j, z - k);
	Vector v3 = Vector(x - (i + 1), y - (j + 1), z - k);
	Vector v4 = Vector(x - i, y - (j + 1), z - k);
	Vector v5 = Vector(x - i, y - j, z - (k + 1));
	Vector v6 = Vector(x - (i + 1), y - j, z - (k + 1));
	Vector v7 = Vector(x - (i + 1), y - (j + 1), z - (k + 1));
	Vector v8 = Vector(x - i, y - (j + 1), z - (k + 1));

	float value = (edge1.dot(v1) * fade(v1.x) * fade(v1.y) * fade(v1.z)) +
		(edge2.dot(v2) * fade(v2.x) * fade(v2.y) * fade(v2.z)) +
		(edge3.dot(v3) * fade(v3.x) * fade(v3.y) * fade(v3.z)) +
		(edge4.dot(v4) * fade(v4.x) * fade(v4.y) * fade(v4.z)) +
		(edge5.dot(v5) * fade(v5.x) * fade(v5.y) * fade(v5.z)) +
		(edge6.dot(v6) * fade(v6.x) * fade(v6.y) * fade(v6.z)) +
		(edge7.dot(v7) * fade(v7.x) * fade(v7.y) * fade(v7.z)) +
		(edge8.dot(v8) * fade(v8.x) * fade(v8.y) * fade(v8.z));


	if (apperance == Appearance::Patch)
	{
		value = (value + 1.0f) / 2.0f;
	}
	else
	{
		value = abs(value);
	}

	return Vector(value, value, value);
}

Vector Texture::get_color(float u, float v)
{
	bool repeat = false;
	if (apperance == Appearance::Repeat)
	{
		u = u - floor(u);
		v = v - floor(v);
		repeat = true;
	}
	else if (apperance == Appearance::Clamp)
	{
		u = std::max(std::min(u, 1.0f), 0.0f);
		v = std::max(std::min(v, 1.0f), 0.0f);
	}

	float i = u * width;
	float j = v * height;

	Vector color = Vector(0.0f);

	if (interpolation == "bilinear")
	{
		color = bilinear_interpolation(i, j, repeat);
	}
	else if (interpolation == "nearest")
	{
		color =  nearest_neighbor(i, j);
	}

	if (degamma)
	{
		color.x = std::pow(color.x / 255, 2.2) * 255;
		color.y = std::pow(color.y / 255, 2.2) * 255;
		color.z = std::pow(color.z / 255, 2.2) * 255;
	}

	return color;
}

Vector Texture::get_perlin_bump(float x, float y, float z)
{
	float epsilon = 0.001f;

	float d = get_perlin_color(x, y, z).x;

	float dx = (get_perlin_color(x + epsilon, y, z).x - d) / epsilon;
	float dy = (get_perlin_color(x, y + epsilon, z).y - d) / epsilon;
	float dz = (get_perlin_color(x, y, z + epsilon).z - d) / epsilon;

	return Vector(dx, dy, dz);
}

Vector Texture::get_bump(float u, float v)
{
	if (apperance == Appearance::Repeat)
	{
		u = u - floor(u);
		v = v - floor(v);
	}
	else if (apperance == Appearance::Clamp)
	{
		u = std::max(std::min(u, 1.0f), 0.0f);
		v = std::max(std::min(v, 1.0f), 0.0f);
	}

	float i = u * width;
	float j = v * height;

	int x = floor(i);
	int y = floor(j);

	int index_a = (y * width + x) * number_of_components;
	int index_b = (y * width + ((x + 1) % width)) * number_of_components;
	int index_c = (((y + 1) % height) * width + x) * number_of_components;

	float a = (image_data[index_a] + image_data[index_a + 1] + image_data[index_a + 2]) / 3.0f;
	float b = (image_data[index_b] + image_data[index_b + 1] + image_data[index_b + 2]) / 3.0f;
	float c = (image_data[index_c] + image_data[index_c + 1] + image_data[index_c + 2]) / 3.0f;

	return Vector(b - a, c - a, 0.0f) * bump_map_multiplier;
}

void Texture::set_decal_mode(std::string mode)
{
	if (mode == "replace_kd")
	{
		decal_mode = DecalMode::ReplaceKD;
	}
	else if (mode == "blend_kd")
	{
		decal_mode = DecalMode::BlendKD;
	}
	else if (mode == "replace_all")
	{
		decal_mode = DecalMode::ReplaceAll;
	}
	else
	{
		decal_mode = DecalMode::None;
	}
}

void Texture::set_appearance(std::string mode)
{
	if (mode == "repeat")
	{
		apperance = Appearance::Repeat;
	}
	else if (mode == "clamp")
	{
		apperance = Appearance::Clamp;
	}
	else if (mode == "patch")
	{
		apperance = Appearance::Patch;
	}
	else if (mode == "vein")
	{
		apperance = Appearance::Vein;
	}
}