#include "Angle.h"

#include "Energy.h"
#include "Gradient.h"

namespace classical {

	Angle::Angle(int atom1, int atom2, int atom3, double degrees, double equilibriumDegrees, double springConstant) 
		: atom1(atom1), atom2(atom2), atom3(atom3), 
		degrees(degrees), equilibriumDegrees(equilibriumDegrees), springConstant(springConstant) {

	}

	void Angle::CalculateEnergy() {
		energy = GetEAngle(degrees, equilibriumDegrees, springConstant);
	}

	void Angle::CalculateGradientMagnitude() {
		gradientMagnitude = GetGMagnitudeAngle(degrees, equilibriumDegrees, springConstant);
	}

	std::ostream& operator<<(std::ostream &stream, const Angle &angle) {
		stream << "Angle: [" << std::endl;
		stream << "\tAtom 1 index: " << angle.atom1 << std::endl;
		stream << "\tAtom 2 index: " << angle.atom2 << std::endl;
		stream << "\tAtom 3 index: " << angle.atom3 << std::endl;
		stream << "\tDegrees: " << angle.degrees << std::endl;
		stream << "\tEquilibrium degrees: " << angle.equilibriumDegrees << std::endl;
		stream << "\tSpring constant: " << angle.springConstant << std::endl;
		stream << "\tEnergy: " << angle.energy << std::endl;
		stream << "\tGradient magnitude: " << angle.gradientMagnitude << std::endl;
		stream << "]";
		return stream;
	}

}