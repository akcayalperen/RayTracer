#include "rtVector.h"
#include <cmath>

Vector::Vector()
{
	x = 0;
	y = 0;
	z = 0;
	w = 0;
}

Vector::Vector(float a)
{
	x = a;
	y = a;
	z = a;
	w = 0;
}

Vector::Vector(float a, float b, float c)
{
	x = a;
	y = b;
	z = c;
	w = 0;
}

Vector::Vector(float a, float b, float c, float d)
{
	x = a;
	y = b;
	z = c;
	w = d;
}

float Vector::length()
{
	return std::sqrt(x * x + y * y + z * z);
}

float Vector::length_squared()
{
	return x * x + y * y + z * z;
}

Vector Vector::cross(const Vector &other)
{
	return Vector(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}

float Vector::dot(const Vector &other)
{
	return x * other.x + y * other.y + z * other.z;
}

Vector Vector::normalized()
{
	return *this / length();
}

void Vector::clamp(float min, float max)
{
	x = std::fminf(std::fmaxf(x, min), max);
	y = std::fminf(std::fmaxf(y, min), max);
	z = std::fminf(std::fmaxf(z, min), max);
}

Vector Vector::find_orthonormal_v()
{
	Vector vec;

	if (y == 0 && x == 0)
	{
		vec = Vector(1, 0, 0);
	}
	else if (y == 0 && z == 0)
	{
		vec = Vector(0, 1, 0);
	}
	else if (x == 0 && z == 0)
	{
		vec = Vector(1, 0, 0);
	}
	else if (x < y && x < z)
	{
		vec = Vector(1, -z, y);
	}
	else if (y < x && y < z)
	{
		vec = Vector(-z, 1, x);
	}
	else
	{
		vec = Vector(-y, x, 1);
	}

	vec = vec.normalized();
	return cross(vec).normalized();
	/*if (fabs(x) > 0.1f)
	{
		vec = Vector(0, 1, 0);
	}
	else
	{
		vec = Vector(1, 0, 0);
	}

	return vec.cross(w).normalized();*/
}

float Vector::determinant_2(float a, float b, float c, float d)
{
	return a * d - b * c;
}

float Vector::determinant_3(const Vector& column1, const Vector& column2, const Vector& column3)
{
	return column1.x * determinant_2(column2.y, column3.y, column2.z, column3.z) -
		column2.x * determinant_2(column1.y, column3.y, column1.z, column3.z) +
		column3.x * determinant_2(column1.y, column2.y, column1.z, column2.z);
}