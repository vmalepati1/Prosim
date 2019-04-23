#include "Torsion.h"

#include "Energy.h"
#include "Gradient.h"

namespace classical {

	Torsion::Torsion(int atom1, int atom2, int atom3, int atom4, 
		double degrees, double halfBarrierHeight, double barrierOffset, 
		int barrierFrequency, int paths) 

		: atom1(atom1), atom2(atom2), atom3(atom3), atom4(atom4), 
		degrees(degrees), halfBarrierHeight(halfBarrierHeight), barrierOffset(barrierOffset), 
		barrierFrequency(barrierFrequency), paths(paths) {

	}

	void Torsion::CalculateEnergy() {
		energy = GetETorsion(degrees, halfBarrierHeight, barrierOffset, barrierFrequency, paths);
	}

	void Torsion::CalculateGradientMagnitude() {
		gradientMagnitude = GetGMagnitudeTorsion(degrees, halfBarrierHeight, barrierOffset, barrierFrequency, paths);
	}

	std::ostream& operator<<(std::ostream &stream, const Torsion &torsion) {
		stream << "Torsion: [" << std::endl;
		stream << "\tAtom 1 index: " << torsion.atom1 << std::endl;
		stream << "\tAtom 2 index: " << torsion.atom2 << std::endl;
		stream << "\tAtom 3 index: " << torsion.atom3 << std::endl;
		stream << "\tAtom 4 index: " << torsion.atom4 << std::endl;
		stream << "\tDegrees: " << torsion.degrees << std::endl;
		stream << "\tRotation half barrier height: " << torsion.halfBarrierHeight << std::endl;
		stream << "\tBarrier offset: " << torsion.barrierOffset << std::endl;
		stream << "\tBarrier frequency: " << torsion.barrierFrequency << std::endl;
		stream << "\tPaths: " << torsion.paths << std::endl;
		stream << "\tEnergy: " << torsion.energy << std::endl;
		stream << "\tGradient magnitude: " << torsion.gradientMagnitude << std::endl;
		stream << "]";
		return stream;
	}

}