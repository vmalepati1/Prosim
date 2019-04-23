#pragma once

#include "Molecule.h"

namespace classical {

	class PQRMolecule : public Molecule {
	public:
		PQRMolecule(const String &pqrFilePath, ForceField *forceField, bool additionalTopologyCalculation);
		~PQRMolecule();

		void CalculateEnergy(const String &kineticType = "none") override;
		void CalculateGradient(const String &gradientType = "analytic") override;
		void CalculateAnalyticGradient() override;
		void CalculateNumericalGradient() override;
		void UpdateInternals() override;
		void CalculateTemperature() override;
		void CalculatePressure() override;
		void CalculateVolume() override;
	private:
		void ReadInPQR();
		void ResolvePQRTokens(const std::vector<String> &tokens);

		void CalculateGNumerical();
	private:
		String m_pqrFilePath;
		ForceField *m_forceField;
	};

}