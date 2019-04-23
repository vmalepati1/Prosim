#include "OutOfPlane.h"

#include "Energy.h"
#include "Gradient.h"

namespace classical {

	OutOfPlane::OutOfPlane(int atom1, int atom2, int atom3, int atom4, double degrees, double halfBarrierHeight)
		: atom1(atom1), atom2(atom2), atom3(atom3), atom4(atom4), degrees(degrees), halfBarrierHeight(halfBarrierHeight) {

	}

	void OutOfPlane::CalculateEnergy() {
		energy = GetEOutOfPlane(degrees, halfBarrierHeight);
	}

	void OutOfPlane::CalculateGradientMagnitude() {
		gradientMagnitude = GetGMagnitudeOutOfPlane(degrees, halfBarrierHeight);
	}

	std::ostream& operator<<(std::ostream &stream, const OutOfPlane &outOfPlane) {
		stream << "OutOfPlane: [" << std::endl;
		stream << "\tAtom 1 index: " << outOfPlane.atom1 << std::endl;
		stream << "\tAtom 2 index: " << outOfPlane.atom2 << std::endl;
		stream << "\tAtom 3 index: " << outOfPlane.atom3 << std::endl;
		stream << "\tAtom 4 index: " << outOfPlane.atom4 << std::endl;
		stream << "\tDegrees: " << outOfPlane.degrees << std::endl;
		stream << "\tRotation half barrier height: " << outOfPlane.halfBarrierHeight << std::endl;
		stream << "\tEnergy: " << outOfPlane.energy << std::endl;
		stream << "\tGradient magnitude: " << outOfPlane.gradientMagnitude << std::endl;
		stream << "]";
		return stream;
	}

}