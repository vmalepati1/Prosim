#pragma once

#include <iostream>
#include <math.h>

namespace classical {

	namespace math {

		struct Vec2;
		struct Vec4;

#pragma pack(push,1)
		struct Vec3 {
			float x, y, z;

			Vec3();
			Vec3(float scalar);
			Vec3(float x, float y, float z);
			Vec3(const Vec2 &other);
			Vec3(float x, float y);
			Vec3(const Vec4 &other);

			Vec3& Add(const Vec3 &other);
			Vec3& Subtract(const Vec3 &other);
			Vec3& Multiply(const Vec3 &other);
			Vec3& Divide(const Vec3 &other);

			Vec3& Add(float other);
			Vec3& Subtract(float other);
			Vec3& Multiply(float other);
			Vec3& Divide(float other);

			friend Vec3& operator+(Vec3 &left, const Vec3 &right);
			friend Vec3& operator-(Vec3 &left, const Vec3 &right);
			friend Vec3& operator*(Vec3 &left, const Vec3 &right);
			friend Vec3& operator/(Vec3 &left, const Vec3 &right);

			friend Vec3 operator+(Vec3 left, float right);
			friend Vec3 operator-(Vec3 left, float right);
			friend Vec3 operator*(Vec3 left, float right);
			friend Vec3 operator/(Vec3 left, float right);

			bool operator==(const Vec3& other);
			bool operator!=(const Vec3& other);

			Vec3& operator+=(const Vec3 &other);
			Vec3& operator-=(const Vec3 &other);
			Vec3& operator*=(const Vec3 &other);
			Vec3& operator/=(const Vec3 &other);

			bool operator<(const Vec3& other) const;
			bool operator<=(const Vec3& other) const;
			bool operator>(const Vec3& other) const;
			bool operator>=(const Vec3& other) const;

			float &operator[](int i);
			const float &operator[](int i) const;

			float Distance(const Vec3 &other) const;

			friend std::ostream& operator<<(std::ostream &stream, const Vec3 &vector);

			float Dot(const Vec3& other) const;
			Vec3 Cross(const Vec3& other) const;

			Vec3 Normalize() const;
		};

	}
#pragma pack(pop)

}

