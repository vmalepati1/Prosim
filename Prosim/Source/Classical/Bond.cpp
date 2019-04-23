#include "Bond.h"

#include "Energy.h"
#include "Gradient.h"

namespace classical {

	Bond::Bond(int atom1, int atom2, double distance, double equilibriumDistance, double springConstant) 
		: atom1(atom1), atom2(atom2), distance(distance), 
		equilibriumDistance(equilibriumDistance), springConstant(springConstant) {

	}

	void Bond::CalculateEnergy() {
		energy = GetEBond(distance, equilibriumDistance, springConstant);
	}

	void Bond::CalculateGradientMagnitude() {
		gradientMagnitude = GetGMagnitudeBond(distance, equilibriumDistance, springConstant);
	}

	std::ostream& operator<<(std::ostream &stream, const Bond &bond) {
		stream << "Bond: [" << std::endl;
		stream << "\tAtom 1 index: " << bond.atom1 << std::endl;
		stream << "\tAtom 2 index: " << bond.atom2 << std::endl;
		stream << "\tDistance: " << bond.distance << std::endl;
		stream << "\tEquilibrium distance: " << bond.equilibriumDistance << std::endl;
		stream << "\tSpring constant: " << bond.springConstant << std::endl;
		stream << "\tEnergy: " << bond.energy << std::endl;
		stream << "\tGradient magnitude: " << bond.gradientMagnitude << std::endl;
		stream << "]";
		return stream;
	}

}