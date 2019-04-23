#pragma once

#include <fstream>

#include "Molecule.h"
#include "SimulationParameters.h"

namespace classical {

	class Simulation {
	public:
		Simulation(Molecule *molecule, const SimulationParameters &simulationParameters);
	protected:
		virtual void OpenOutputFiles() = 0;
		virtual void CloseOutputFiles();
		virtual void FlushBuffers();
		virtual void WriteGeometry() = 0;
		virtual void WriteValue(int totalFloatChars, int decimalChars, double value, char printType = 'f', int nSpace = 1);
		virtual void WriteEnergyTerms(int totalFloatChars, int decimalChars, char printType) = 0;
		virtual void WriteEnergy() = 0;
		virtual void WriteEnergyHeader() = 0;
		virtual void PrintStatus() = 0;
	protected:
		Molecule *m_molecule;
		SimulationParameters m_parameters;
		std::ofstream m_energyFile;
		std::ofstream m_geometryFile;
	};

}
