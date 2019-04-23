#include "PQRMolecule.h"

#include <algorithm>
#include <fstream>

#include "Constants.h"
#include "Energy.h"
#include "Geometry.h"
#include "Gradient.h"
#include "Topology.h"

namespace classical {

	PQRMolecule::PQRMolecule(const String &pqrFilePath, ForceField *forceField, bool additionalTopologyCalculation)
		: m_pqrFilePath(pqrFilePath), m_forceField(forceField) {

		ReadInPQR();

		CalculateBondGraphFromBonds(m_nAtoms, m_bondGraph, m_bonds);

		if (additionalTopologyCalculation) {
			CalculateBondGraph(m_atoms, m_bondGraph, m_forceField);
		}

		std::cout << "Calculated bond graph" << std::endl;

		CalculateBonds(m_atoms, m_bondGraph, m_bonds, m_forceField);
		std::cout << "Calculated bonds" << std::endl;

		CalculateAngles(m_atoms, m_bondGraph, m_angles, m_forceField);
		std::cout << "Calculated angles" << std::endl;

		CalculateTorsions(m_atoms, m_bondGraph, m_torsions, m_forceField);
		std::cout << "Calculated torsions" << std::endl;

		CalculateOutOfPlanes(m_atoms, m_bondGraph, m_outOfPlanes, m_forceField);
		std::cout << "Calculated out-of-planes" << std::endl;

		CalculateNonInts(m_bonds, m_angles, m_torsions, m_nonInts);
		std::cout << "Calculated non-interacting pairs" << std::endl;

		m_nBonds = m_bonds.size();
		m_nAngles = m_angles.size();
		m_nTorsions = m_torsions.size();
		m_nOutOfPlanes = m_outOfPlanes.size();

		m_dielectric = 1.0;
		m_mass = 0.0;
		m_kBox = 250.0;
		m_boundary = 1.0E10;
		m_boundaryType = "sphere";
		m_origin = math::Vec3();
		m_volume = INFINITY;
		m_temperature = 0.0;
		m_pressure = 0.0;
		m_virial = 0.0;

		for (int i = 0; i < m_nAtoms; i++) {
			m_gBonds.push_back(math::Vec3());
			m_gAngles.push_back(math::Vec3());
			m_gTorsions.push_back(math::Vec3());
			m_gOutOfPlanes.push_back(math::Vec3());
			m_gVDW.push_back(math::Vec3());
			m_gElst.push_back(math::Vec3());
			m_gBound.push_back(math::Vec3());
			m_gBonded.push_back(math::Vec3());
			m_gNonBonded.push_back(math::Vec3());
			m_gPotential.push_back(math::Vec3());
			m_gKinetic.push_back(math::Vec3());
			m_gTotal.push_back(math::Vec3());
		}
	}

	PQRMolecule::~PQRMolecule() {
		for (int i = 0; i < m_nAtoms; i++) {
			delete m_atoms[i];
		}

		for (int i = 0; i < m_nBonds; i++) {
			delete m_bonds[i];
		}

		for (int i = 0; i < m_nAngles; i++) {
			delete m_angles[i];
		}

		for (int i = 0; i < m_nTorsions; i++) {
			delete m_torsions[i];
		}

		for (int i = 0; i < m_nOutOfPlanes; i++) {
			delete m_outOfPlanes[i];
		}
	}

	void PQRMolecule::CalculateEnergy(const String &kineticType) {
		m_eBonds = GetEBonds(m_bonds);
		m_eAngles = GetEAngles(m_angles);
		m_eTorsions = GetETorsions(m_torsions);
		m_eOutOfPlanes = GetEOutOfPlanes(m_outOfPlanes);

#ifdef PS_OPTIMIZED
		math::Vec2 nonBondedEnergy = m_nonBondedGPUCalculator.GetENonBonded(m_atoms, m_nonInts, m_dielectric);
#else
		math::Vec2 nonBondedEnergy = GetENonBonded(m_atoms, m_nonInts, m_dielectric);
#endif

		m_eVDW = nonBondedEnergy.x;
		m_eElst = nonBondedEnergy.y;

		m_eBound = GetEBound(m_atoms, m_kBox, m_boundary, m_origin, m_boundaryType);
		m_eKinetic = GetEKinetic(m_atoms, kineticType);
		m_eBonded = m_eBonds + m_eAngles + m_eTorsions + m_eOutOfPlanes;

		m_eNonBonded = m_eVDW + m_eElst;

		m_ePotential = m_eBonded + m_eNonBonded + m_eBound;

		m_eTotal = m_ePotential + m_eKinetic;
	}

	void PQRMolecule::CalculateGradient(const String &gradientType) {
		if (gradientType == "analytic") {
			CalculateAnalyticGradient();
		}
		else if (gradientType == "numerical") {
			CalculateNumericalGradient();
		} else {
			std::cout << "Unknown gradient type: " << gradientType << std::endl;
			std::cout << "Use 'analytic' or 'numerical'" << std::endl;
			return;
		}

		std::fill(m_gBonded.begin(), m_gBonded.end(), 0);
		std::fill(m_gNonBonded.begin(), m_gNonBonded.end(), 0);
		std::fill(m_gTotal.begin(), m_gTotal.end(), 0);

		for (int i = 0; i < m_nAtoms; i++) {
			m_gBonded[i].Add(m_gBonds[i]);
			m_gBonded[i].Add(m_gAngles[i]);
			m_gBonded[i].Add(m_gTorsions[i]);
			m_gBonded[i].Add(m_gOutOfPlanes[i]);

			m_gNonBonded[i].Add(m_gVDW[i]);
			m_gNonBonded[i].Add(m_gElst[i]);

			m_gTotal[i].Add(m_gBonded[i]);
			m_gTotal[i].Add(m_gNonBonded[i]);
			m_gTotal[i].Add(m_gBound[i]);
		}
	}

	void PQRMolecule::CalculateAnalyticGradient() {
		CalculateGBonds(m_gBonds, m_bonds, m_atoms);
		CalculateGAngles(m_gAngles, m_angles, m_atoms, m_bondGraph);
		CalculateGTorsions(m_gTorsions, m_torsions, m_atoms, m_bondGraph);
		CalculateGOutOfPlanes(m_gOutOfPlanes, m_outOfPlanes, m_atoms, m_bondGraph);
	}

	void PQRMolecule::CalculateNumericalGradient() {
		CalculateGNumerical();
	}

	void PQRMolecule::UpdateInternals() {
		UpdateBonds(m_bonds, m_atoms, m_bondGraph);
		UpdateAngles(m_angles, m_atoms, m_bondGraph);
		UpdateTorsions(m_torsions, m_atoms, m_bondGraph);
		UpdateOutOfPlanes(m_outOfPlanes, m_atoms, m_bondGraph);
	}

	void PQRMolecule::CalculateTemperature() {
		m_temperature = GetTemperature(m_eKinetic, m_nAtoms);
	}

	void PQRMolecule::CalculatePressure() {
		m_virial = GetVirial(m_gTotal, m_atoms);
		m_pressure = GetPressure(m_atoms, m_temperature, m_virial, m_volume);
	}

	void PQRMolecule::CalculateVolume() {
		m_volume = GetVolume(m_boundary, m_boundaryType);
	}

	void PQRMolecule::ReadInPQR() {
		std::ifstream file(m_pqrFilePath);

		if (file.is_open()) {
			String line;

			while (std::getline(file, line)) {
				if (line.find("REMARK") != 0 && !line.empty()) {
					std::vector<String> tokens = utils::Tokenize(line);

					ResolvePQRTokens(tokens);
				}
			}

			file.close();
		}
		else {
			std::cout << "Could not open PQR file " << m_pqrFilePath << std::endl;
			return;
		}

		m_nAtoms = m_atoms.size();
	}

	void PQRMolecule::ResolvePQRTokens(const std::vector<String> &tokens) {

		if (tokens[0] == "ATOM" || tokens[0] == "HETATM") {
			if (tokens.size() == 10 || tokens.size() == 11) {
				bool isShortFormat = (tokens.size() == 11 ? false : true);

				int serial = utils::NextInt(tokens[1]);
				String atomName = tokens[2];
				String residueName = tokens[3];
				int chainID = (isShortFormat ? 0 : utils::NextInt(tokens[4]));
				int residueNumber = (isShortFormat ? utils::NextInt(tokens[4]) : utils::NextInt(tokens[5]));
				double x = (isShortFormat ? utils::ToDouble(tokens[5]) : utils::ToDouble(tokens[6]));
				double y = (isShortFormat ? utils::ToDouble(tokens[6]) : utils::ToDouble(tokens[7]));
				double z = (isShortFormat ? utils::ToDouble(tokens[7]) : utils::ToDouble(tokens[8]));
				double charge = (isShortFormat ? utils::ToDouble(tokens[8]) : utils::ToDouble(tokens[9]));
				double radius = (isShortFormat ? utils::ToDouble(tokens[9]) : utils::ToDouble(tokens[10]));

				math::Vec3 position = math::Vec3(x, y, z);

				m_atoms.push_back(new Atom(atomName, position, charge, m_forceField));
			}
		}

		if (tokens[0] == "CONECT") {
			switch (tokens.size() - 1) {
			case 2: {
				/* Normal two-atom bond. */

				int atom1Index = utils::NextInt(tokens[1]) - 1;
				int atom2Index = utils::NextInt(tokens[2]) - 1;

				Atom *atom1 = m_atoms[atom1Index];
				Atom *atom2 = m_atoms[atom2Index];

				double distance = GetRij(atom1->position, atom2->position);
				double equilibriumDistance = m_forceField->GetBondEquilibriumLength(atom1->type, atom2->type);
				double springConstant = m_forceField->GetBondSpringConstant(atom1->type, atom2->type);

				m_bonds.push_back(new Bond(atom1Index, atom2Index, distance, equilibriumDistance, springConstant));

				break;
			}
			default:
				break;
			}
		}
	}

	void PQRMolecule::CalculateGNumerical() {
		std::fill(m_gBonds.begin(), m_gBonds.end(), 0);
		std::fill(m_gAngles.begin(), m_gAngles.end(), 0);
		std::fill(m_gTorsions.begin(), m_gTorsions.end(), 0);
		std::fill(m_gOutOfPlanes.begin(), m_gOutOfPlanes.end(), 0);
		std::fill(m_gVDW.begin(), m_gVDW.end(), 0);
		std::fill(m_gElst.begin(), m_gElst.end(), 0);
		std::fill(m_gBound.begin(), m_gBound.end(), 0);

		for (int i = 0; i < m_nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				double q = m_atoms[i]->position[j];

				double qp = q + 0.5 * NUMERICAL_DISPLACEMENT;

				m_atoms[i]->position[j] = qp;

				UpdateInternals();
				CalculateEnergy();

				double epBond = m_eBonds;
				double epAngle = m_eAngles;
				double epTorsion = m_eTorsions;
				double epOutOfPlane = m_eOutOfPlanes;
				double epVDW = m_eVDW;
				double epElst = m_eElst;
				double epBound = m_eBound;

				double qm = q - 0.5 * NUMERICAL_DISPLACEMENT;
				
				UpdateInternals();
				CalculateEnergy();

				double emBond = m_eBonds;
				double emAngle = m_eAngles;
				double emTorsion = m_eTorsions;
				double emOutOfPlane = m_eOutOfPlanes;
				double emVDW = m_eVDW;
				double emElst = m_eElst;
				double emBound = m_eBound;

				double displacement = qp - qm;

				m_atoms[i]->position[j] = q;
				m_gBonds[i][j] = (epBond - emBond) / displacement;
				m_gAngles[i][j] = (epAngle - emAngle) / displacement;
				m_gTorsions[i][j] = (epTorsion - emTorsion) / displacement;
				m_gOutOfPlanes[i][j] = (epOutOfPlane - emOutOfPlane) / displacement;
				m_gVDW[i][j] = (epVDW - emVDW) / displacement;
				m_gElst[i][j] = (epElst - emElst) / displacement;
				m_gBound[i][j] = (epBound - emBound) / displacement;
			}
		}

		UpdateInternals();
	}

}