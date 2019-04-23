#pragma once

#include <sstream>

namespace classical {

	struct Bond {
		Bond(int atom1, int atom2, double distance, double equilibriumDistance, double springConstant);

		void CalculateEnergy();
		void CalculateGradientMagnitude();

		friend std::ostream& operator<<(std::ostream &stream, const Bond &bond);

		int atom1;
		int atom2;
		double distance;
		double equilibriumDistance;
		double springConstant;
		double energy;
		double gradientMagnitude;
	};

}