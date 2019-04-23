#pragma once

#include <sstream>

namespace classical {

	struct Torsion {
		Torsion(int atom1, int atom2, int atom3, int atom4, double degrees, double halfBarrierHeight, double barrierOffset, int barrierFrequency, int paths);

		void CalculateEnergy();
		void CalculateGradientMagnitude();

		friend std::ostream& operator<<(std::ostream &stream, const Torsion &torsion);

		int atom1;
		int atom2;
		int atom3;
		int atom4;
		double degrees;
		double halfBarrierHeight;
		double barrierOffset;
		int barrierFrequency;
		int paths;
		double energy;
		double gradientMagnitude;
	};

}