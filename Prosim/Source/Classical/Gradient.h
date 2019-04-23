#pragma once

#include "Angle.h"
#include "Atom.h"
#include "Bond.h"
#include "OutOfPlane.h"
#include "Torsion.h"

#include "Math/PSMath.h"
#include "Utils/String.h"

namespace classical {

	double GetGMagnitudeBond(double rij, double req, double kb);
	double GetGMagnitudeAngle(double aijk, double aeq, double ka);
	double GetGMagnitudeTorsion(double tijkl, double vn, double gamma, double nfold, double paths);
	double GetGMagnitudeOutOfPlane(double oijkl, double vn);
	double GetGMagnitudeVDWIJ(double rij, double epsij, double roij);
	double GetGMagnitudeElstIJ(double rij, double qi, double qj, double epsilon);
	math::Vec3 GetGMagnitudeBoundI(double kBox, double bound, const math::Vec3 &position, const math::Vec3 &origin, const String &boundType);
	std::pair<math::Vec3, math::Vec3> GetGDirectionInteraction(const math::Vec3 &position1, const math::Vec3 &position2, double r12 = -1);
	std::tuple<math::Vec3, math::Vec3, math::Vec3> GetGDirectionAngle(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, double r21 = -1, double r23 = -1);
	std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> GetGDirectionTorsion(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, const math::Vec3 &position4, double r12 = -1, double r23 = -1, double r34 = -1);
	std::tuple<math::Vec3, math::Vec3, math::Vec3, math::Vec3> GetGDirectionOutOfPlane(const math::Vec3 &position1, const math::Vec3 &position2, const math::Vec3 &position3, const math::Vec3 &position4, double degrees, double r31 = -1, double r32 = -1, double r34 = -1);
	void CalculateGBonds(std::vector<math::Vec3> &gBonds, const std::vector<Bond *> &bonds, const std::vector<Atom *> &atoms);
	void CalculateGAngles(std::vector<math::Vec3> &gAngles, const std::vector<Angle *> &angles, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph);
	void CalculateGTorsions(std::vector<math::Vec3> &gTorsions, const std::vector<Torsion *> &torsions, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph);
	void CalculateGOutOfPlanes(std::vector<math::Vec3> &gOutOfPlanes, const std::vector<OutOfPlane *> &outOfPlanes, const std::vector<Atom *> atoms, std::map<int, std::map<int, double>> &bondGraph);
	void CalculateGNonBonded(std::vector<math::Vec3> &gVDW, std::vector<math::Vec3> &gElst, const std::vector<Atom *> &atoms, const std::vector<std::pair<int, int>> &nonInts, double dielectric);
	void CalculateGBound(std::vector<math::Vec3> &gBound, const std::vector<Atom *> &atoms, double kBox, double bound, const math::Vec3 &origin, const String &boundType);
	double GetVirial(std::vector<math::Vec3> &gTotal, const std::vector<Atom *> &atoms);
	double GetPressure(const std::vector<Atom *> &atoms, double temperature, double virial, double volume);

}