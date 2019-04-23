#pragma once

#include <iostream>
#include <math.h>

#include "Vec3.h"

namespace classical {

	namespace math {

		struct Vec2 {
			float x, y;

			Vec2();
			Vec2(float scalar);
			Vec2(float x, float y);
			Vec2(const Vec3 &vector);

			Vec2& Add(const Vec2 &other);
			Vec2& Subtract(const Vec2 &other);
			Vec2& Multiply(const Vec2 &other);
			Vec2& Divide(const Vec2 &other);

			friend Vec2& operator+(Vec2 &left, const Vec2 &right);
			friend Vec2& operator-(Vec2 &left, const Vec2 &right);
			friend Vec2& operator*(Vec2 &left, const Vec2 &right);
			friend Vec2& operator/(Vec2 &left, const Vec2 &right);

			bool operator==(const Vec2& other);
			bool operator!=(const Vec2& other);

			Vec2& operator+=(const Vec2 &other);
			Vec2& operator-=(const Vec2 &other);
			Vec2& operator*=(const Vec2 &other);
			Vec2& operator/=(const Vec2 &other);

			bool operator<(const Vec2& other) const;
			bool operator<=(const Vec2& other) const;
			bool operator>(const Vec2& other) const;
			bool operator>=(const Vec2& other) const;

			float &operator[](int i);
			const float &operator[](int i) const;

			float Magnitude() const;
			Vec2 Normalize() const;
			float Distance(const Vec2 &other) const;
			float Dot(const Vec2 &other) const;

			friend std::ostream& operator<<(std::ostream &stream, const Vec2 &vector);

		};

	}

}
