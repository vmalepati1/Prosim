#pragma once

#include <CL/cl.h>

#include "Angle.h"
#include "Atom.h"
#include "Bond.h"
#include "OutOfPlane.h"
#include "Torsion.h"

#include "Math/PSMath.h"
#include "Utils/String.h"

namespace classical {

	class NonBondedEnergyGPUCalculator {
	public:
		NonBondedEnergyGPUCalculator();
		~NonBondedEnergyGPUCalculator();

		math::Vec2 GetENonBonded(const std::vector<Atom *> &atoms, std::vector<int> &nonInts, double dielectric);
	private:
		cl_kernel kernel;
		cl_program program;
		cl_mem kCombinationMatrixMem;
		cl_mem kNonIntsMem;
		cl_mem kNonIntsSizeMem;
		cl_mem kPositionsMem;
		cl_mem kChargesMem;
		cl_mem kVdwRadiiMem;
		cl_mem kVdwAttractionMagnitudesMem;
		cl_mem kDielectricMem;
		cl_mem kCeuToKCalMem;
		cl_mem keVDWMem;
		cl_mem keElstMem;
		cl_command_queue commandQueue;
		cl_context context;
	};

	double GetEBond(double rij, double req, double kb);
	double GetEAngle(double aijk, double aeq, double ka);
	double GetETorsion(double tijkl, double vn, double gamma, int nfold, int paths);
	double GetEOutOfPlane(double oijkl, double vn);
	double GetEVDWIJ(double rij, double epsij, double roij);
	double GetEElstIJ(double rij, double qi, double qj, double epsilon);
	double GetEBoundI(double kBox, double bound, const math::Vec3 &position, const math::Vec3 &origin, const String &boundType);
	double GetEKineticI(double mass, const math::Vec3 &velocity);
	double GetEBonds(const std::vector<Bond *> &bonds);
	double GetEAngles(const std::vector<Angle *> &angles);
	double GetETorsions(const std::vector<Torsion *> &torsions);
	double GetEOutOfPlanes(const std::vector<OutOfPlane *> &outOfPlanes);
	math::Vec2 GetENonBonded(const std::vector<Atom *> &atoms, std::vector<int> &nonInts, double dielectric);
	double GetEBound(const std::vector<Atom *> &atoms, double kBox, double boundary, const math::Vec3 &origin, const String &boundType);
	double GetEKinetic(const std::vector<Atom *> &atoms, const String &kineticType = "none");
	double GetTemperature(double eKinetic, int natoms);

}