#pragma once

#include <sstream>

namespace classical {

	struct Angle {
		Angle(int atom1, int atom2, int atom3, double degrees, double equilibriumDegrees, double springConstant);

		void CalculateEnergy();
		void CalculateGradientMagnitude();

		friend std::ostream& operator<<(std::ostream &stream, const Angle &angle);

		int atom1;
		int atom2;
		int atom3;
		double degrees;
		double equilibriumDegrees;
		double springConstant;
		double energy;
		double gradientMagnitude;
	};

}