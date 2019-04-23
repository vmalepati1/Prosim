#include "Atom.h"

#include <algorithm>

namespace classical {

	Atom::Atom(const String &type, const math::Vec3 &position, double charge, ForceField *forceField, double ro, double epsilon) 
		: type(type), position(position), charge(charge), vdwRadius(ro), vdwAttractionMagnitude(epsilon) {

		InferElementFromType();

		mass = forceField->GetAtomicMass(element);
		covalentRadius = forceField->GetCovalentRadius(element);

		velocity = math::Vec3();
		acceleration = math::Vec3();
		previousVelocity = math::Vec3();
		previousAcceleration = math::Vec3();
	}

	Atom::Atom(const String &type, const math::Vec3 &position, double charge, ForceField *forceField) 
		: type(type), position(position), charge(charge) {

		vdwRadius = forceField->GetVanDerWaalsRadius(type);
		vdwAttractionMagnitude = forceField->GetVanDerWaalsAttractionMagnitude(type);

		InferElementFromType();

		mass = forceField->GetAtomicMass(element);
		covalentRadius = forceField->GetCovalentRadius(element);

		velocity = math::Vec3();
		acceleration = math::Vec3();
		previousVelocity = math::Vec3();
		previousAcceleration = math::Vec3();
	}

	void Atom::InferElementFromType() {
		if (type.size() == 1 || !islower(type[1])) {
			element = type.substr(0, 1);
		}
		else {
			element = type.substr(0, 2);
		}

		std::transform(element.begin(), element.end(), element.begin(), ::toupper);
	}

	std::ostream& operator<<(std::ostream &stream, const Atom &atom) {
		stream << "Atom: [" << std::endl;
		stream << "\tPosition: " << atom.position << std::endl;
		stream << "\tVelocity: " << atom.velocity << std::endl;
		stream << "\tAcceleration: " << atom.acceleration << std::endl;
		stream << "\tPrevious velocity: " << atom.velocity << std::endl;
		stream << "\tPrevious acceleration: " << atom.acceleration << std::endl;
		stream << "\tElement: " << atom.element << std::endl;
		stream << "\tType: " << atom.type << std::endl;
		stream << "\tCharge: " << atom.charge << std::endl;
		stream << "\tVan der Waals radius: " << atom.vdwRadius << std::endl;
		stream << "\tVan der Waals epsilon: " << atom.vdwAttractionMagnitude << std::endl;
		stream << "\tAtomic mass: " << atom.mass << std::endl;
		stream << "\tCovalent radius: " << atom.covalentRadius << std::endl;
		stream << "]";
		return stream;
	}

}