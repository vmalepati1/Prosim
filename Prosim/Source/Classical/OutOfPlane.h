#pragma once

#include <sstream>

namespace classical {

	struct OutOfPlane {
		OutOfPlane(int atom1, int atom2, int atom3, int atom4, double degrees, double halfBarrierHeight);
	
		void CalculateEnergy();
		void CalculateGradientMagnitude();

		friend std::ostream& operator<<(std::ostream &stream, const OutOfPlane &outOfPlane);

		int atom1;
		int atom2;
		int atom3;
		int atom4;
		double degrees;
		double halfBarrierHeight;
		double energy;
		double gradientMagnitude;
	};

}