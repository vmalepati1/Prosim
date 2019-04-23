#include "Vec3.h"
#include "Vec2.h"
#include "Vec4.h"

namespace classical {

	namespace math {

		Vec3::Vec3() {
			this->x = 0.0f;
			this->y = 0.0f;
			this->z = 0.0f;
		}

		Vec3::Vec3(float scalar) {
			this->x = scalar;
			this->y = scalar;
			this->z = scalar;
		}

		Vec3::Vec3(float x, float y, float z) {
			this->x = x;
			this->y = y;
			this->z = z;
		}

		Vec3::Vec3(const Vec2 &other) {
			this->x = other.x;
			this->y = other.y;
			this->z = 0.0f;
		}

		Vec3::Vec3(float x, float y) {
			this->x = x;
			this->y = y;
			this->z = 0.0f;
		}

		Vec3::Vec3(const Vec4 &other) {
			this->x = other.x;
			this->y = other.y;
			this->z = other.z;
		}

		Vec3& Vec3::Add(const Vec3 &other) {
			x += other.x;
			y += other.y;
			z += other.z;
			return *this;
		}

		Vec3& Vec3::Subtract(const Vec3 &other) {
			x -= other.x;
			y -= other.y;
			z -= other.z;
			return *this;
		}

		Vec3& Vec3::Multiply(const Vec3 &other) {
			x *= other.x;
			y *= other.y;
			z *= other.z;
			return *this;
		}

		Vec3& Vec3::Divide(const Vec3 &other) {
			x /= other.x;
			y /= other.y;
			z /= other.z;
			return *this;
		}
		Vec3& Vec3::Add(float other)
		{
			x += other;
			y += other;
			z += other;

			return *this;
		}

		Vec3& Vec3::Subtract(float other)
		{
			x -= other;
			y -= other;
			z -= other;

			return *this;
		}

		Vec3& Vec3::Multiply(float other)
		{
			x *= other;
			y *= other;
			z *= other;

			return *this;
		}

		Vec3& Vec3::Divide(float other)
		{
			x /= other;
			y /= other;
			z /= other;

			return *this;
		}


		Vec3& operator+(Vec3 &left, const Vec3 &right) {
			return left.Add(right);
		}

		Vec3& operator-(Vec3 &left, const Vec3 &right) {
			return left.Subtract(right);
		}

		Vec3& operator*(Vec3 &left, const Vec3 &right) {
			return left.Multiply(right);
		}

		Vec3& operator/(Vec3 &left, const Vec3 &right) {
			return left.Divide(right);
		}

		Vec3 operator+(Vec3 left, float right) {
			return left.Add(right);
		}

		Vec3 operator-(Vec3 left, float right) {
			return left.Subtract(right);
		}

		Vec3 operator*(Vec3 left, float right) {
			return left.Multiply(right);
		}

		Vec3 operator/(Vec3 left, float right) {
			return left.Divide(right);
		}

		Vec3& Vec3::operator+=(const Vec3 &other) {
			*this = *this + other;
			return *this;
		}

		Vec3& Vec3::operator-=(const Vec3 &other) {
			*this = *this - other;
			return *this;
		}

		Vec3& Vec3::operator*=(const Vec3 &other) {
			*this = *this * other;
			return *this;
		}

		Vec3& Vec3::operator/=(const Vec3 &other) {
			*this = *this / other;
			return *this;
		}

		bool Vec3::operator<(const Vec3& other) const
		{
			return x < other.x && y < other.y && z < other.z;
		}

		bool Vec3::operator<=(const Vec3& other) const
		{
			return x <= other.x && y <= other.y && z <= other.z;
		}

		bool Vec3::operator>(const Vec3& other) const
		{
			return x > other.x && y > other.y && z > other.z;
		}

		bool Vec3::operator>=(const Vec3& other) const
		{
			return x >= other.x && y >= other.y && z >= other.z;
		}

		float &Vec3::operator[](int i) {
			if (i == 0) return x;
			else if (i == 1) return y;
			else if (i == 2) return z;
			else {
				std::cout << "Vec3 index " << i << " out of range (returned x)!" << std::endl;
				return x;
			}
		}

		const float &Vec3::operator[](int i) const {
			if (i == 0) return x;
			else if (i == 1) return y;
			else if (i == 2) return z;
			else {
				std::cout << "Vec3 index " << i << " out of range (returned x)!" << std::endl;
				return -1;
			}
		}

		bool Vec3::operator==(const Vec3& other) {
			return x == other.x && y == other.y && z == other.z;
		}

		bool Vec3::operator!=(const Vec3& other) {
			return !(*this == other);
		}

		float Vec3::Distance(const Vec3 &other) const {
			float a = x - other.x;
			float b = y - other.y;
			float c = z - other.z;
			return sqrt(a * a + b * b + c * c);
		}

		std::ostream& operator<<(std::ostream &stream, const Vec3 &vector) {
			stream << "Vec3 (" << vector.x << ", " << vector.y << ", " << vector.z << ")";
			return stream;
		}

		float Vec3::Dot(const Vec3& other) const
		{
			return x * other.x + y * other.y + z * other.z;
		}

		Vec3 Vec3::Cross(const Vec3& other) const
		{
			return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
		}

		Vec3 Vec3::Normalize() const {
			float magnitude = sqrt((x * x) + (y * y) + (z * z));
			return Vec3(x / magnitude, y / magnitude, z / magnitude);
		}

	}

}