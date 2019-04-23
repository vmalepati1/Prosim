#pragma once

#include <string>
#include <vector>

#include "Angle.h"
#include "Atom.h"
#include "Bond.h"
#include "Energy.h"
#include "ForceField.h"
#include "OutOfPlane.h"
#include "Torsion.h"

namespace classical {

	class Simulation;
	class MolecularDynamics;

	class Molecule {
	public:
		virtual void CalculateEnergy(const String &kineticType = "none") = 0;
		virtual void CalculateGradient(const String &gradientType = "analytic") = 0;
		virtual void CalculateAnalyticGradient() = 0;
		virtual void CalculateNumericalGradient() = 0;
		virtual void UpdateInternals() = 0;
		virtual void CalculateTemperature() = 0;
		virtual void CalculatePressure() = 0;
		virtual void CalculateVolume() = 0;

		/* Some getters/setters contain 'Member' in them as there are conflicting method names in the global function space. */

		inline double GetKBox() const { return m_kBox; }
		inline double GetBoundary() const { return m_boundary; }
		inline const String &GetBoundaryType() const { return m_boundaryType; }
		inline const math::Vec3 &GetOrigin() const { return m_origin; }

		inline void SetKBox(double kBox) { m_kBox = kBox; }
		inline void SetBoundary(double boundary) { m_boundary = boundary; }
		inline void SetBoundaryType(const String &boundaryType) { m_boundaryType = boundaryType; }
		inline void SetOrigin(const math::Vec3 &origin) { m_origin = origin; }

		inline double GetDielectric() const { return m_dielectric; }
		inline double GetMass() const { return m_mass; }
		inline double GetMemberVolume() const { return m_volume; }
		inline double GetMemberTemperature() const { return m_temperature; }
		inline double GetMemberPressure() const { return m_pressure; }
		inline double GetMemberVirial() const { return m_virial; }

		inline double SetDielectric(double dielectric) { m_dielectric = dielectric; }

		inline const std::vector<Atom *> &GetAtoms() const { return m_atoms; }
		inline const std::vector<Bond *> &GetBonds() const { return m_bonds; }
		inline const std::vector<Angle *> &GetAngles() const { return m_angles; }
		inline const std::vector<Torsion *> &GetTorsions() const { return m_torsions; }
		inline const std::vector<OutOfPlane *> &GetOutOfPlanes() const { return m_outOfPlanes; }

		inline const std::vector<int> &GetNonInts() const { return m_nonInts; }
		inline const std::map<int, std::map<int, double>> &GetBondGraph() const { return m_bondGraph; }

		inline int GetNAtoms() const { return m_nAtoms; }
		inline int GetNBonds() const { return m_nBonds; }
		inline int GetNAngles() const { return m_nAngles; }
		inline int GetNTorsions() const { return m_nTorsions; }
		inline int GetNOutOfPlanes() const { return m_nOutOfPlanes; }

		inline double GetMemberEBonds() const { return m_eBonds; }
		inline double GetMemberEAngles() const { return m_eAngles; }
		inline double GetMemberETorsions() const { return m_eTorsions; }
		inline double GetMemberEOutOfPlanes() const { return m_eOutOfPlanes; }
		inline double GetMemberEVDW() const { return m_eVDW; }
		inline double GetMemberEElst() const { return m_eElst; }
		inline double GetMemberEBound() const { return m_eBound; }
		inline double GetMemberEBonded() const { return m_eBonded; }
		inline double GetMemberENonBonded() const { return m_eNonBonded; }
		inline double GetMemberEPotential() const { return m_ePotential; }
		inline double GetMemberEKinetic() const { return m_eKinetic; }
		inline double GetMemberETotal() const { return m_eTotal; }

		inline const std::vector<math::Vec3> &GetGBonds() const { return m_gBonds; }
		inline const std::vector<math::Vec3> &GetGAngles() const { return m_gAngles; }
		inline const std::vector<math::Vec3> &GetGTorsions() const { return m_gTorsions; }
		inline const std::vector<math::Vec3> &GetGOutOfPlanes() const { return m_gOutOfPlanes; }
		inline const std::vector<math::Vec3> &GetGVDW() const { return m_gVDW; }
		inline const std::vector<math::Vec3> &GetGElst() const { return m_gElst; }
		inline const std::vector<math::Vec3> &GetGBound() const { return m_gBound; }
		inline const std::vector<math::Vec3> &GetGBonded() const { return m_gBonded; }
		inline const std::vector<math::Vec3> &GetGNonBonded() const { return m_gNonBonded; }
		inline const std::vector<math::Vec3> &GetGPotential() const { return m_gPotential; }
		inline const std::vector<math::Vec3> &GetGKinetic() const { return m_gKinetic; }
		inline const std::vector<math::Vec3> &GetGTotal() const { return m_gTotal; }

	protected:
		double m_kBox;
		double m_boundary;
		String m_boundaryType;
		math::Vec3 m_origin;

		double m_dielectric;
		double m_mass;
		double m_volume;
		double m_temperature;
		double m_pressure;
		double m_virial;

		std::vector<Atom *> m_atoms;
		std::vector<Bond *> m_bonds;
		std::vector<Angle *> m_angles;
		std::vector<Torsion *> m_torsions;
		std::vector<OutOfPlane *> m_outOfPlanes;

		std::vector<int> m_nonInts;
		std::map<int, std::map<int, double>> m_bondGraph;

		int m_nAtoms;
		int m_nBonds;
		int m_nAngles;
		int m_nTorsions;
		int m_nOutOfPlanes;

		/* All energies are in kcal/mol. */
		double m_eBonds;
		double m_eAngles;
		double m_eTorsions;
		double m_eOutOfPlanes;
		double m_eVDW;
		double m_eElst;
		double m_eBound;
		double m_eBonded;
		double m_eNonBonded;
		double m_ePotential;
		double m_eKinetic;
		double m_eTotal;

		std::vector<math::Vec3> m_gBonds;
		std::vector<math::Vec3> m_gAngles;
		std::vector<math::Vec3> m_gTorsions;
		std::vector<math::Vec3> m_gOutOfPlanes;
		std::vector<math::Vec3> m_gVDW;
		std::vector<math::Vec3> m_gElst;
		std::vector<math::Vec3> m_gBound;
		std::vector<math::Vec3> m_gBonded;
		std::vector<math::Vec3> m_gNonBonded;
		std::vector<math::Vec3> m_gPotential;
		std::vector<math::Vec3> m_gKinetic;
		std::vector<math::Vec3> m_gTotal;

		NonBondedEnergyGPUCalculator m_nonBondedGPUCalculator;

		friend class Simulation;
		friend class MolecularDynamics;
	};

}