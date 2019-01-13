#pragma once
#ifndef __RT_MATRIX__H_
#define __RT_MATRIX__H_

#include <stdexcept>
#include "rtVector.h"
#include <math.h>

#define PI 3.14159265f

class Matrix
{
public:
	float m[16];

	Matrix()
	{
		m[0] = m[5] = m[10] = m[15] = 1;
		m[1] = m[2] = m[3] = m[4] = 0;
		m[6] = m[7] = m[8] = m[9] = 0;
		m[11] = m[12] = m[13] = m[14] = 0;
	}
	Matrix(const float mat[16])
	{
		m[0] = mat[0]; m[1] = mat[1]; m[2] = mat[2]; m[3] = mat[3];
		m[4] = mat[4]; m[5] = mat[5]; m[6] = mat[6]; m[7] = mat[7];
		m[8] = mat[8]; m[9] = mat[9]; m[10] = mat[10]; m[11] = mat[11];
		m[12] = mat[12]; m[13] = mat[13]; m[14] = mat[14]; m[15] = mat[15];
	}
	Matrix(float m0, float m1, float m2, float m3,
		float m4, float m5, float m6, float m7,
		float m8, float m9, float m10, float m11,
		float m12, float m13, float m14, float m15)
	{
		m[0] = m0; m[1] = m1; m[2] = m2; m[3] = m3;
		m[4] = m4; m[5] = m5; m[6] = m6; m[7] = m7;
		m[8] = m8; m[9] = m9; m[10] = m10; m[11] = m11;
		m[12] = m12; m[13] = m13; m[14] = m14; m[15] = m15;
	}

	static Matrix Identity()
	{
		return Matrix();
	}

	static Matrix TranslationMatrix(Vector vec)
	{
		Matrix matrix = Matrix::Identity();

		matrix.m[3] = vec.x;
		matrix.m[7] = vec.y;
		matrix.m[11] = vec.z;

		return matrix;
	}

	static Matrix ScaleMatrix(Vector vec)
	{
		Matrix matrix = Matrix::Identity();

		matrix.m[0] = vec.x;
		matrix.m[5] = vec.y;
		matrix.m[10] = vec.z;

		return matrix;
	}

	static Matrix RotationMatrix(Vector vec, float angle = 0)
	{
		Matrix matrix = Matrix::Identity();

		angle = angle * PI / 180.0f;
		float cosine = cos(angle);
		float sine = sin(angle);

		Vector u = vec.normalized();

		float x = u.x;
		float y = u.y;
		float z = u.z;
		float t = 1 - cosine;

		matrix.m[0] = t * x * x + cosine;
		matrix.m[1] = t * x * y - sine * z;
		matrix.m[2] = t * x * z + sine * y;
		matrix.m[3] = 0.0;

		matrix.m[4] = t * x * y + sine * z;
		matrix.m[5] = t * y * y + cosine;
		matrix.m[6] = t * y * z - sine * x;
		matrix.m[7] = 0.0;

		matrix.m[8] = t * x * z - sine * y;
		matrix.m[9] = t * y * z + sine * x;
		matrix.m[10] = t * z * z + cosine;
		matrix.m[11] = 0.0;

		matrix.m[12] = 0.0;
		matrix.m[13] = 0.0;
		matrix.m[14] = 0.0;
		matrix.m[15] = 1.0;

		return matrix;
	}

	Matrix operator * (const Matrix& mat) const
	{
		Matrix result;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				result.m[i * 4 + j] = m[i * 4] * mat.m[j] +
					m[i * 4 + 1] * mat.m[j + 4] +
					m[i * 4 + 2] * mat.m[j + 8] +
					m[i * 4 + 3] * mat.m[j + 12];
			}
		}
		return result;
	}

	Vector transform(const Vector& vec) const
	{
		Vector result;
		result.x = m[0] * vec.x + m[1] * vec.y + m[2] * vec.z + m[3] * vec.w;
		result.y = m[4] * vec.x + m[5] * vec.y + m[6] * vec.z + m[7] * vec.w;
		result.z = m[8] * vec.x + m[9] * vec.y + m[10] * vec.z + m[11] * vec.w;
		result.w = m[12] * vec.x + m[13] * vec.y + m[14] * vec.z + m[15] * vec.w;
		return result;
	}

	Vector transform_point(const Vector& vec) const
	{
		Vector result;
		result.x = m[0] * vec.x + m[1] * vec.y + m[2] * vec.z + m[3];
		result.y = m[4] * vec.x + m[5] * vec.y + m[6] * vec.z + m[7];
		result.z = m[8] * vec.x + m[9] * vec.y + m[10] * vec.z + m[11];
		return result;
	}

	Vector transform_vector(const Vector& vec) const
	{
		Vector result;
		result.x = m[0] * vec.x + m[1] * vec.y + m[2] * vec.z;
		result.y = m[4] * vec.x + m[5] * vec.y + m[6] * vec.z;
		result.z = m[8] * vec.x + m[9] * vec.y + m[10] * vec.z;
		return result;
	}

	Matrix inverse()
	{
		float inv[16];
		float det;

		inv[0] = m[5] * m[10] * m[15] -
			m[5] * m[11] * m[14] -
			m[9] * m[6] * m[15] +
			m[9] * m[7] * m[14] +
			m[13] * m[6] * m[11] -
			m[13] * m[7] * m[10];

		inv[4] = -m[4] * m[10] * m[15] +
			m[4] * m[11] * m[14] +
			m[8] * m[6] * m[15] -
			m[8] * m[7] * m[14] -
			m[12] * m[6] * m[11] +
			m[12] * m[7] * m[10];

		inv[8] = m[4] * m[9] * m[15] -
			m[4] * m[11] * m[13] -
			m[8] * m[5] * m[15] +
			m[8] * m[7] * m[13] +
			m[12] * m[5] * m[11] -
			m[12] * m[7] * m[9];

		inv[12] = -m[4] * m[9] * m[14] +
			m[4] * m[10] * m[13] +
			m[8] * m[5] * m[14] -
			m[8] * m[6] * m[13] -
			m[12] * m[5] * m[10] +
			m[12] * m[6] * m[9];

		inv[1] = -m[1] * m[10] * m[15] +
			m[1] * m[11] * m[14] +
			m[9] * m[2] * m[15] -
			m[9] * m[3] * m[14] -
			m[13] * m[2] * m[11] +
			m[13] * m[3] * m[10];

		inv[5] = m[0] * m[10] * m[15] -
			m[0] * m[11] * m[14] -
			m[8] * m[2] * m[15] +
			m[8] * m[3] * m[14] +
			m[12] * m[2] * m[11] -
			m[12] * m[3] * m[10];

		inv[9] = -m[0] * m[9] * m[15] +
			m[0] * m[11] * m[13] +
			m[8] * m[1] * m[15] -
			m[8] * m[3] * m[13] -
			m[12] * m[1] * m[11] +
			m[12] * m[3] * m[9];

		inv[13] = m[0] * m[9] * m[14] -
			m[0] * m[10] * m[13] -
			m[8] * m[1] * m[14] +
			m[8] * m[2] * m[13] +
			m[12] * m[1] * m[10] -
			m[12] * m[2] * m[9];

		inv[2] = m[1] * m[6] * m[15] -
			m[1] * m[7] * m[14] -
			m[5] * m[2] * m[15] +
			m[5] * m[3] * m[14] +
			m[13] * m[2] * m[7] -
			m[13] * m[3] * m[6];

		inv[6] = -m[0] * m[6] * m[15] +
			m[0] * m[7] * m[14] +
			m[4] * m[2] * m[15] -
			m[4] * m[3] * m[14] -
			m[12] * m[2] * m[7] +
			m[12] * m[3] * m[6];

		inv[10] = m[0] * m[5] * m[15] -
			m[0] * m[7] * m[13] -
			m[4] * m[1] * m[15] +
			m[4] * m[3] * m[13] +
			m[12] * m[1] * m[7] -
			m[12] * m[3] * m[5];

		inv[14] = -m[0] * m[5] * m[14] +
			m[0] * m[6] * m[13] +
			m[4] * m[1] * m[14] -
			m[4] * m[2] * m[13] -
			m[12] * m[1] * m[6] +
			m[12] * m[2] * m[5];

		inv[3] = -m[1] * m[6] * m[11] +
			m[1] * m[7] * m[10] +
			m[5] * m[2] * m[11] -
			m[5] * m[3] * m[10] -
			m[9] * m[2] * m[7] +
			m[9] * m[3] * m[6];

		inv[7] = m[0] * m[6] * m[11] -
			m[0] * m[7] * m[10] -
			m[4] * m[2] * m[11] +
			m[4] * m[3] * m[10] +
			m[8] * m[2] * m[7] -
			m[8] * m[3] * m[6];

		inv[11] = -m[0] * m[5] * m[11] +
			m[0] * m[7] * m[9] +
			m[4] * m[1] * m[11] -
			m[4] * m[3] * m[9] -
			m[8] * m[1] * m[7] +
			m[8] * m[3] * m[5];

		inv[15] = m[0] * m[5] * m[10] -
			m[0] * m[6] * m[9] -
			m[4] * m[1] * m[10] +
			m[4] * m[2] * m[9] +
			m[8] * m[1] * m[6] -
			m[8] * m[2] * m[5];

		det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

		if (det == 0)
		{
			throw std::runtime_error("Error: Determinant is zero.");
		}

		det = 1.0f / det;

		for (int i = 0; i < 16; ++i)
		{
			inv[i] *= det;
		}

		return Matrix(inv);
	}

	Matrix transpose()
	{
		float transpose[16];

		transpose[0] = m[0]; transpose[1] = m[4]; transpose[2] = m[8]; transpose[3] = m[12];
		transpose[4] = m[1]; transpose[5] = m[5]; transpose[6] = m[9]; transpose[7] = m[13];
		transpose[8] = m[2]; transpose[9] = m[6]; transpose[10] = m[10]; transpose[11] = m[14];
		transpose[12] = m[3]; transpose[13] = m[7]; transpose[14] = m[11]; transpose[15] = m[15];

		return Matrix(transpose);
	}
};

#endif