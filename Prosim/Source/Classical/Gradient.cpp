#include "Gradient.h"

#include <math.h>

#include "Constants.h"
#include "Geometry.h"

#include "Utils/IterationTools.h"

namespace classical {

	double GetGMagnitudeBond(double rij, double req, double kb) {
		return 2.0 * kb * (rij - req);
	}

	double GetGMagnitudeAngle(double aijk, double aeq, double ka) {
		return 2.0 * ka * (DEGREES_TO_RADIANS * (aijk - aeq));
	}

	double GetGMagnitudeTorsion(double tijkl, double vn, double gamma, double nfold, double paths) {
		return -vn * nfold * sin(DEGREES_TO_RADIANS * (nfold * tijkl - gamma)) / paths;
	}

	double GetGMagnitudeOutOfPlane(double oijkl, double vn) {
		return -2.0 * vn * sin(DEGREES_TO_RADIANS * (2.0 * oijkl - 180.0));
	}

	double GetGMagnitudeVDWIJ(double rij, double epsij, double roij) {
		double rrelij = roij / rij;
		return -12.0 * (epsij / roij) * (pow(rrelij, 13) - pow(rrelij, 7));
	}

	double GetGMagnitudeElstIJ(double rij, double qi, double qj, double epsilon) {
		return -CEU_TO_KCAL * qi * qj / (epsilon * pow(rij, 2));
	}

	math::Vec3 GetGMagnitudeBoundI(double kBox, double bound, const math::Vec3 &position, const math::Vec3 &origin, const String &boundType) {
		math::Vec3 gBoundI;

		if (boundType == "cube") {
			double sign = (position.x - origin.x) <= 0 ? 1.0 : -1.0;
			double scale = abs(position.x - origin.x);
			gBoundI.x = (-2.0 * sign * scale * kBox * (abs(position.x) - bound));

			sign = (position.y - origin.y) <= 0 ? 1.0 : -1.0;
			scale = abs(position.y - origin.y);
			gBoundI.x = (-2.0 * sign * scale * kBox * (abs(position.y) - bound));

			sign = (position.z - origin.z) <= 0 ? 1.0 : -1.0;
			scale = abs(position.z - origin.z);
			gBoundI.x = (-2.0 * sign * scale * kBox * (abs(position.z) - bound));
		}
		else if (boundType == "sphere") {
			double rio = GetRij(position, origin);
			math::Vec3 uio = GetUij(origin, position);
			double scale = (double)rio >= bound;
			gBoundI = uio * 2.0 * scale * kBox * (rio - bound);
		}

		return gBoundI;
	}

	std::pair<math::Vec3, math::Vec3> GetGDirectionInteraction(const math::Vec3 &position1, const math::Vec3 &position2, double r12) {
		math::Vec3 gDir1 = GetUij(position2, position1, r12);
		math::Vec3 gDir2 = math::Vec3(gDir1) * -1.0;
		return std::make_pair(gDir1, gDir2);
	}

	std::tuple<math::Vec3, math::Vec3, math::Vec3> GetGDirectionAngle(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, double r21, double r23) {
		math::Vec3 u21 = GetUij(position2, position1, r21);
		math::Vec3 u23 = GetUij(position2, position3, r23);
		math::Vec3 cp = GetUcp(u21, u23);
		math::Vec3 gDir1 = GetUcp(u21, cp) / r21;
		math::Vec3 gDir3 = GetUcp(cp, u23) / r23;
		math::Vec3 gDir2 = (math::Vec3(gDir1) + math::Vec3(gDir3)) * -1.0;
		return std::make_tuple(gDir1, gDir2, gDir3);
	}

	std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> GetGDirectionTorsion(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, const math::Vec3 &position4, double r12, double r23, double r34) {
		math::Vec3 u21 = GetUij(position2, position1, r12);
		math::Vec3 u34 = GetUij(position3, position4, r34);
		math::Vec3 u23 = GetUij(position2, position3, r23);
		math::Vec3 u32 = math::Vec3(u23) * -1.0;
		double a123 = GetAijk(position1, position2, position3, r12, r23);
		double a432 = GetAijk(position4, position3, position2, r34, r23);
		double s123 = sin(DEGREES_TO_RADIANS * a123);
		double s432 = sin(DEGREES_TO_RADIANS * a432);
		double c123 = cos(DEGREES_TO_RADIANS * a123);
		double c432 = cos(DEGREES_TO_RADIANS * a432);
		math::Vec3 gDir1 = GetUcp(u21, u23) / (r12 * s123);
		math::Vec3 gDir4 = GetUcp(u34, u32) / (r34 * s432);
		math::Vec3 gDir2 = math::Vec3(gDir1) * (r12 / r23 * c123 - 1.0) - math::Vec3(gDir4) * (r34 / r23 * c432);
		math::Vec3 gDir3 = math::Vec3(gDir4) * (r34 / r23 * c432 - 1.0) - math::Vec3(gDir1) * (r12 / r23 * c123);
		return std::make_tuple(gDir1, gDir2, gDir3, gDir4);
	}

	std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> GetGDirectionOutOfPlane(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, const math::Vec3 &position4, double degrees, double r31, double r32, double r34) {
		math::Vec3 u31 = GetUij(position3, position1, r31);
		math::Vec3 u32 = GetUij(position3, position2, r32);
		math::Vec3 u34 = GetUij(position3, position4, r34);
		math::Vec3 cp3234 = GetCp(u32, u34);
		math::Vec3 cp3431 = GetCp(u34, u31);
		math::Vec3 cp3132 = GetCp(u31, u32);
		double a132 = GetAijk(position1, position3, position2);
		double s132 = sin(DEGREES_TO_RADIANS * a132);
		double c132 = cos(DEGREES_TO_RADIANS * a132);
		double cOOP = cos(DEGREES_TO_RADIANS * degrees);
		double tOOP = tan(DEGREES_TO_RADIANS * degrees);
		math::Vec3 gDir1 = ((math::Vec3(u31) - math::Vec3(u32) * c132)) * (1.0 / r31) * (math::Vec3(cp3234) / (cOOP * s132) - (pow(tOOP / s132, 2)));
		math::Vec3 gDir2 = ((math::Vec3(cp3431) / (cOOP * s132) - (pow(tOOP / s132, 2)) - (math::Vec3(u32) - math::Vec3(u31) * c132)) * (1.0 / r32));
		math::Vec3 gDir4 = ((cp3132 / (cOOP * s132) - (u34 * tOOP)) * (1.0 / r34));
		math::Vec3 gDir3 = (math::Vec3(gDir1) + math::Vec3(gDir2) + math::Vec3(gDir4)) * -1.0;
		return std::make_tuple(gDir1, gDir2, gDir3, gDir4);
	}

	void CalculateGBonds(std::vector<math::Vec3> &gBonds, const std::vector<Bond *> &bonds, const std::vector<Atom *> &atoms) {
		std::fill(gBonds.begin(), gBonds.end(), 0);

		for (Bond *bond : bonds) {
			bond->CalculateGradientMagnitude();
			math::Vec3 p1 = atoms[bond->atom1]->position;
			math::Vec3 p2 = atoms[bond->atom2]->position;
			std::tuple<math::Vec3, math::Vec3> directions = GetGDirectionInteraction(p1, p2, bond->distance);
			math::Vec3 dir1 = std::get<0>(directions);
			math::Vec3 dir2 = std::get<1>(directions);
			gBonds[bond->atom1] += dir1 * bond->gradientMagnitude;
			gBonds[bond->atom2] += dir2 * bond->gradientMagnitude;
		}
	}

	void CalculateGAngles(std::vector<math::Vec3> &gAngles, const std::vector<Angle *> &angles, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph) {
		std::fill(gAngles.begin(), gAngles.end(), 0);

		for (Angle *angle : angles) {
			angle->CalculateGradientMagnitude();
			math::Vec3 p1 = atoms[angle->atom1]->position;
			math::Vec3 p2 = atoms[angle->atom2]->position;
			math::Vec3 p3 = atoms[angle->atom3]->position;
			double r12 = bondGraph[angle->atom1][angle->atom2];
			double r23 = bondGraph[angle->atom2][angle->atom3];
			std::tuple<math::Vec3, math::Vec3, math::Vec3> directions = GetGDirectionAngle(p1, p2, p3, r12, r23);
			math::Vec3 dir1 = std::get<0>(directions);
			math::Vec3 dir2 = std::get<1>(directions);
			math::Vec3 dir3 = std::get<2>(directions);
			gAngles[angle->atom1] += dir1 * angle->gradientMagnitude;
			gAngles[angle->atom2] += dir2 * angle->gradientMagnitude;
			gAngles[angle->atom3] += dir3 * angle->gradientMagnitude;
		}
	}

	void CalculateGTorsions(std::vector<math::Vec3> &gTorsions, const std::vector<Torsion *> &torsions, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph) {
		std::fill(gTorsions.begin(), gTorsions.end(), 0);

		for (Torsion *torsion : torsions) {
			torsion->CalculateGradientMagnitude();
			math::Vec3 p1 = atoms[torsion->atom1]->position;
			math::Vec3 p2 = atoms[torsion->atom2]->position;
			math::Vec3 p3 = atoms[torsion->atom3]->position;
			math::Vec3 p4 = atoms[torsion->atom4]->position;
			double r12 = bondGraph[torsion->atom1][torsion->atom2];
			double r23 = bondGraph[torsion->atom2][torsion->atom3];
			double r34 = bondGraph[torsion->atom3][torsion->atom4];
			std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> directions = GetGDirectionTorsion(p1, p2, p3, p4, r12, r23, r34);
			math::Vec3 dir1 = std::get<0>(directions);
			math::Vec3 dir2 = std::get<1>(directions);
			math::Vec3 dir3 = std::get<2>(directions);
			math::Vec3 dir4 = std::get<3>(directions);
			gTorsions[torsion->atom1] += dir1 * torsion->gradientMagnitude;
			gTorsions[torsion->atom2] += dir2 * torsion->gradientMagnitude;
			gTorsions[torsion->atom3] += dir3 * torsion->gradientMagnitude;
			gTorsions[torsion->atom4] += dir4 * torsion->gradientMagnitude;
		}
	}

	void CalculateGOutOfPlanes(std::vector<math::Vec3> &gOutOfPlanes, const std::vector<OutOfPlane *> &outOfPlanes, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph) {
		std::fill(gOutOfPlanes.begin(), gOutOfPlanes.end(), 0);

		for (OutOfPlane *outOfPlane : outOfPlanes) {
			outOfPlane->CalculateGradientMagnitude();
			math::Vec3 p1 = atoms[outOfPlane->atom1]->position;
			math::Vec3 p2 = atoms[outOfPlane->atom2]->position;
			math::Vec3 p3 = atoms[outOfPlane->atom3]->position;
			math::Vec3 p4 = atoms[outOfPlane->atom4]->position;
			double r31 = bondGraph[outOfPlane->atom3][outOfPlane->atom1];
			double r32 = bondGraph[outOfPlane->atom3][outOfPlane->atom2];
			double r34 = bondGraph[outOfPlane->atom3][outOfPlane->atom4];
			std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> directions = GetGDirectionOutOfPlane(p1, p2, p3, p4, outOfPlane->degrees, r31, r32, r34);
			math::Vec3 dir1 = std::get<0>(directions);
			math::Vec3 dir2 = std::get<1>(directions);
			math::Vec3 dir3 = std::get<2>(directions);
			math::Vec3 dir4 = std::get<3>(directions);
			gOutOfPlanes[outOfPlane->atom1] += dir1 * outOfPlane->gradientMagnitude;
			gOutOfPlanes[outOfPlane->atom2] += dir2 * outOfPlane->gradientMagnitude;
			gOutOfPlanes[outOfPlane->atom3] += dir3 * outOfPlane->gradientMagnitude;
			gOutOfPlanes[outOfPlane->atom4] += dir4 * outOfPlane->gradientMagnitude;
		}

	}

	void CalculateGNonBonded(std::vector<math::Vec3> &gVDW, std::vector<math::Vec3> &gElst, const std::vector<Atom *> &atoms, const std::vector<std::pair<int, int>> &nonInts, double dielectric) {
		std::fill(gVDW.begin(), gVDW.end(), 0);
		std::fill(gElst.begin(), gElst.end(), 0);

		int natoms = atoms.size();

		utils::IterationMatrix matrix(utils::CombinationsNR(natoms, 2), 2);

		utils::CombinationKN(matrix, 2, natoms);

#ifdef PS_OPTIMIZED
#pragma omp parallel for
#endif
		for (int n = 0; n < utils::CombinationsNR(natoms, 2); n++) {
			int i = matrix(n, 0);
			int j = matrix(n, 1);

			if (std::find(nonInts.begin(), nonInts.end(), std::make_pair(i, j)) != nonInts.end()) continue;

			Atom *atom1 = atoms[i];
			Atom *atom2 = atoms[j];

			double distance = GetRij(atom1->position, atom2->position);
			std::tuple<math::Vec3, math::Vec3> directions = GetGDirectionInteraction(atom1->position, atom2->position, distance);
			math::Vec3 dir1 = std::get<0>(directions);
			math::Vec3 dir2 = std::get<1>(directions);

			double vdwAttractionMagnitudeIJ = sqrt(atom1->vdwAttractionMagnitude) * sqrt(atom2->vdwAttractionMagnitude);
			double vdwRadiusIJ = atom1->vdwRadius + atom2->vdwRadius;

			double gElstMagnitude = GetGMagnitudeElstIJ(distance, atom1->charge, atom2->charge, dielectric);
			double gVDWMagnitude = GetGMagnitudeVDWIJ(distance, vdwAttractionMagnitudeIJ, vdwRadiusIJ);

			gVDW[i] += math::Vec3(dir1) * gVDWMagnitude;
			gVDW[j] += math::Vec3(dir2) * gVDWMagnitude;
			gElst[i] += math::Vec3(dir1) * gElstMagnitude;
			gElst[j] += math::Vec3(dir2) * gElstMagnitude;
		}
	}

	void CalculateGBound(std::vector<math::Vec3> &gBound, const std::vector<Atom *> &atoms, double kBox, double bound, const math::Vec3 &origin, const String &boundType) {
		std::fill(gBound.begin(), gBound.end(), 0);

		int i = 0;
		for (Atom *atom : atoms) {
			gBound[i] += GetGMagnitudeBoundI(kBox, bound, atom->position, origin, boundType);
			i++;
		}
	}

	double GetVirial(std::vector<math::Vec3> &gTotal, const std::vector<Atom *> &atoms) {
		double virial = 0.0;

		int i = 0;
		for (Atom *atom : atoms) {
			virial += gTotal[i].x;
			virial += gTotal[i].y;
			virial += gTotal[i].z;
			i++;
		}

		return virial;
	}

	double GetPressure(const std::vector<Atom *> &atoms, double temperature, double virial, double volume) {
		return KCAL_A_MOL_TO_PA * (atoms.size() * BOLTZMANN_CONSTANT * temperature + virial / 3) / volume;
	}

}