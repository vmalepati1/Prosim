#include "Vec2.h"

namespace classical {

	namespace math {

		Vec2::Vec2() {
			this->x = 0.0f;
			this->y = 0.0f;
		}

		Vec2::Vec2(float x, float y) {
			this->x = x;
			this->y = y;
		}

		Vec2::Vec2(float scalar) {
			this->x = scalar;
			this->y = scalar;
		}

		Vec2::Vec2(const Vec3 &vector) {
			this->x = vector.x;
			this->y = vector.y;
		}

		Vec2& Vec2::Add(const Vec2 &other) {
			x += other.x;
			y += other.y;
			return *this;
		}

		Vec2& Vec2::Subtract(const Vec2 &other) {
			x -= other.x;
			y -= other.y;
			return *this;
		}

		Vec2& Vec2::Multiply(const Vec2 &other) {
			x *= other.x;
			y *= other.y;
			return *this;
		}

		Vec2& Vec2::Divide(const Vec2 &other) {
			x /= other.x;
			y /= other.y;
			return *this;
		}

		Vec2& operator+(Vec2 &left, const Vec2 &right) {
			return left.Add(right);
		}

		Vec2& operator-(Vec2 &left, const Vec2 &right) {
			return left.Subtract(right);
		}

		Vec2& operator*(Vec2 &left, const Vec2 &right) {
			return left.Multiply(right);
		}

		Vec2& operator/(Vec2 &left, const Vec2 &right) {
			return left.Divide(right);
		}

		Vec2& Vec2::operator+=(const Vec2 &other) {
			*this = *this + other;
			return *this;
		}

		Vec2& Vec2::operator-=(const Vec2 &other) {
			*this = *this - other;
			return *this;
		}

		Vec2& Vec2::operator*=(const Vec2 &other) {
			*this = *this * other;
			return *this;
		}

		Vec2& Vec2::operator/=(const Vec2 &other) {
			*this = *this / other;
			return *this;
		}

		bool Vec2::operator==(const Vec2& other) {
			return x == other.x && y == other.y;
		}

		bool Vec2::operator!=(const Vec2& other) {
			return !(*this == other);
		}

		bool Vec2::operator<(const Vec2& other) const
		{
			return x < other.x && y < other.y;
		}

		bool Vec2::operator<=(const Vec2& other) const
		{
			return x <= other.x && y <= other.y;
		}

		bool Vec2::operator>(const Vec2& other) const
		{
			return x > other.x && y > other.y;
		}

		bool Vec2::operator>=(const Vec2& other) const
		{
			return x >= other.x && y >= other.y;
		}

		float &Vec2::operator[](int i) {
			if (i == 0) return x;
			else if (i == 1) return y;
			else {
				std::cout << "Vec2 index " << i << " out of range (returned x)!" << std::endl;
				return x;
			}
		}

		const float &Vec2::operator[](int i) const {
			if (i == 0) return x;
			else if (i == 1) return y;
			else {
				std::cout << "Vec2 index " << i << " out of range (returned x)!" << std::endl;
				return x;
			}
		}

		float Vec2::Magnitude() const {
			return sqrt(x * x + y * y);
		}

		Vec2 Vec2::Normalize() const {
			float length = Magnitude();
			return Vec2(x / length, y / length);
		}

		float Vec2::Distance(const Vec2 &other) const {
			float a = x - other.x;
			float b = y - other.y;
			return sqrt(a * a + b * b);
		}

		float Vec2::Dot(const Vec2 &other) const {
			return x * other.x + y * other.y;
		}

		std::ostream& operator<<(std::ostream &stream, const Vec2 &vector) {
			stream << "Vec2 (" << vector.x << ", " << vector.y << ")";
			return stream;
		}

	}

}