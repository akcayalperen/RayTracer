#pragma once
#ifndef __RT_VECTOR__H_
#define __RT_VECTOR__H_

class Vector
{
public:
	Vector();
	Vector(float a);
	Vector(float a, float b, float c);
	Vector(float a, float b, float c, float d);

	Vector operator - () const
	{
		return Vector(-x, -y, -z);
	}

	float length();
	float length_squared();
	Vector cross(const Vector &other);
	float dot(const Vector &other);
	Vector normalized();
	void clamp(float min, float max);

	Vector find_orthonormal_v();
	static float determinant_2(float a, float b, float c, float d);
	static float determinant_3(const Vector& a, const Vector& b, const Vector& c);

	union
	{
		float coordinate[4];
		struct
		{
			float x;
			float y;
			float z;
			float w;
		};
	};
};

inline Vector operator+(const Vector& a, const Vector& b)
{
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator-(const Vector& a, const Vector& b)
{
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector operator*(const Vector& a, const float& b)
{
	return Vector(a.x*b, a.y*b, a.z*b);
}

inline Vector operator*(const float& b, const Vector& a)
{
	return Vector(a.x*b, a.y*b, a.z*b);
}

inline Vector operator*(const Vector& a, const Vector& b)
{
	return Vector(a.x*b.x, a.y*b.y, a.z*b.z);
}

inline Vector operator/(const Vector& a, const float& b)
{
	return Vector(a.x / b, a.y / b, a.z / b);
}

#endif