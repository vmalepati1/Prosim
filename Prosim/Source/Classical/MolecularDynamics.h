#pragma once

#include "Simulation.h"

namespace classical {

	class MolecularDynamics : public Simulation {
	public:
		MolecularDynamics(Molecule *molecule, const SimulationParameters &parameters);

		void Run();
	private:
		void OpenOutputFiles() override;
		void WriteGeometry() override;
		void WriteEnergyTerms(int totalFloatChars, int decimalChars, char printType) override;
		void WriteEnergy() override;
		void WriteEnergyHeader() override;
		void PrintStatus() override;

		void InitializeVelocities();
		void EquilibrateTemperature();
		void UpdateAccelerations();
		void UpdateVelocities(double deltaTime);
		void UpdatePositions(double deltaTime);
		void CheckPrint(double timeStep, bool printAll = false);

	private:
		double m_lastTime;
		double m_currentTime;
		double m_eTemperature;

		double m_eTime;
		double m_gTime;
	};

}