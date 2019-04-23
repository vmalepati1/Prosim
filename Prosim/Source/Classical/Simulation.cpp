#include "Simulation.h"

namespace classical {

	Simulation::Simulation(Molecule *molecule, const SimulationParameters &simulationParameters)
		: m_molecule(molecule), m_parameters(simulationParameters) {

		m_molecule->m_kBox = m_parameters.GetBoundarySpring();
		m_molecule->m_boundary = m_parameters.GetBoundary();
		m_molecule->CalculateVolume();
		m_molecule->m_boundaryType = m_parameters.GetBoundaryType();
		m_molecule->CalculateVolume();
		m_molecule->m_origin = m_parameters.GetOrigin();
	}

	void Simulation::CloseOutputFiles() {
		PrintStatus();

		m_energyFile.close();
		m_geometryFile.close();
	}

	void Simulation::FlushBuffers() {
		m_energyFile.flush();
		m_geometryFile.flush();

		fflush(stdout);
	}

	void Simulation::WriteValue(int totalFloatChars, int decimalChars, double value, char printType, int nSpace) {
		if (printType == 'f') {
			m_energyFile << utils::StringWithFormat("%*s%*.*f", nSpace, "", totalFloatChars, decimalChars, value);
		}
		else if (printType == 'e') {
			m_energyFile << utils::StringWithFormat("%*s%*.*e", nSpace, "", totalFloatChars, decimalChars, value);
		}
	}

}