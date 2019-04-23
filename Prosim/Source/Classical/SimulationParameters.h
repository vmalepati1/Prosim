#pragma once

#include "Math/PSMath.h"
#include "Utils/String.h"

namespace classical {

	class Simulation;
	class MolecularDynamics;

	class SimulationParameters {
	public:
		SimulationParameters();
		SimulationParameters(const String &filePath);

		inline const String &GetFilePath() { return m_filePath;  }
		inline double GetDesiredTemperature() { return m_desiredTemperature; }
		inline double GetDesiredPressure() { return m_desiredPressure; }
		inline double GetBoundarySpring() { return m_boundarySpring; }
		inline double GetBoundary() { return m_boundary; }
		inline const String &GetBoundaryType() { return m_boundaryType; }
		inline const math::Vec3 &GetOrigin() { return m_origin; }
		inline double GetTotalTime() { return m_totalTime; }
		inline int GetTotalConfigurations() { return m_totalConfigurations; }
		inline double GetTimeStep() { return m_timeStep; }
		inline double GetGeometryWaitTime() { return m_geometryWaitTime; }
		inline int GetGeometryConfigurations() { return m_geometryConfigurations; }
		inline const String &GetEnergyOutputFilePath() const { return m_energyOutputFilePath; }
		inline const String &GetGeometryOutputFilePath() const { return m_geometryOutputFilePath; }
		inline double GetEnergyWaitTime() { return m_energyWaitTime; }
		inline int GetEnergyConfigurations() { return m_energyConfigurations; }
		inline double GetStatusWaitTime() { return m_statusWaitTime; }
		inline double GetEquilibriumTime() { return m_equilibriumTime; }
		inline double GetEquilibriumRate() { return m_equilibriumRate; }
		inline int GetRandomSeed() { return m_randomSeed; }
		inline int GetEnergyWriteDigits() { return m_energyWriteDigits; }
		inline int GetEnergyWriteChars() { return m_energyWriteChars; }
		inline int GetGeometryWriteDigits() { return m_geometryWriteDigits; }
		inline int GetGeometryWriteChars() { return m_geometryWriteChars; }
		inline int GetTimeWriteDigits() { return m_timeWriteDigits; }
		inline int GetTimeWriteChars() { return m_timeWriteChars; }

		friend std::ostream& operator<<(std::ostream &stream, const SimulationParameters &simulationParameters);
	private:
		void SetDefaultParameters();
		void ResolveParameter(const String &key, const String &value);
	private:
		String m_filePath;

		double m_desiredTemperature;
		double m_desiredPressure;
		double m_boundarySpring;
		double m_boundary;
		String m_boundaryType;
		math::Vec3 m_origin;
		double m_totalTime;
		int m_totalConfigurations;
		double m_timeStep;
		double m_geometryWaitTime;
		int m_geometryConfigurations;
		String m_energyOutputFilePath;
		String m_geometryOutputFilePath;
		double m_energyWaitTime;
		int m_energyConfigurations;
		double m_statusWaitTime;
		double m_equilibriumTime;
		double m_equilibriumRate;
		int m_randomSeed;
		int m_energyWriteDigits;
		int m_energyWriteChars;
		int m_geometryWriteDigits;
		int m_geometryWriteChars;
		int m_timeWriteDigits;
		int m_timeWriteChars;

		friend class Simulation;
		friend class MolecularDynamics;
	};

}