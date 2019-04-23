#pragma once

#include "ForceField.h"
#include "Math/PSMath.h"
#include "Utils/String.h"

namespace classical {

	struct Atom {
		Atom(const String &type, const math::Vec3 &position, double charge, ForceField *forceField, double ro, double epsilon);
		Atom(const String &type, const math::Vec3 &position, double charge, ForceField *forceField);

		friend std::ostream& operator<<(std::ostream &stream, const Atom &atom);

		math::Vec3 position;
		math::Vec3 velocity;
		math::Vec3 acceleration;
		math::Vec3 previousVelocity;
		math::Vec3 previousAcceleration;

		String element;
		String type;

		double charge;
		double vdwRadius;
		double vdwAttractionMagnitude;
		double mass;
		double covalentRadius;

	private:
		void InferElementFromType();
	};

}