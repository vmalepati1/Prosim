#include "SimulationParameters.h"

#include <algorithm>
#include <fstream>
#include <iostream>

namespace classical {

	SimulationParameters::SimulationParameters()
		: m_filePath("") {
		SetDefaultParameters();
	}

	SimulationParameters::SimulationParameters(const String &filePath) 
		: m_filePath(filePath) {
		SetDefaultParameters();

		std::ifstream file(filePath);

		if (file.is_open()) {
			String line;

			while (std::getline(file, line)) {
				if (line.find("#") != 0 && !line.empty() && line.find("=")) {
					String part = line.substr(0, line.find("#"));
					String key = part.substr(0, line.find("="));
					String value = part.substr(line.find("=") + 1);

					ResolveParameter(key, value);
				}
			}

			file.close();
		}
		else {
			std::cout << "Could not open parameter file " << filePath << std::endl;
			return;
		}
	}

	std::ostream& operator<<(std::ostream &stream, const SimulationParameters &simulationParameters) {
		stream << "Parameters: [" << std::endl;
		stream << "\tDesired temperature: " << simulationParameters.m_desiredTemperature << std::endl;
		stream << "\tDesired pressure: " << simulationParameters.m_desiredPressure << std::endl;
		stream << "\tBoundary spring: " << simulationParameters.m_boundarySpring << std::endl;
		stream << "\tBoundary: " << simulationParameters.m_boundary << std::endl;
		stream << "\tBoundary type: " << simulationParameters.m_boundaryType << std::endl;
		stream << "\tOrigin: " << simulationParameters.m_origin << std::endl;
		stream << "\tTotal time: " << simulationParameters.m_totalTime << std::endl;
		stream << "\tTotal configurations: " << simulationParameters.m_totalConfigurations << std::endl;
		stream << "\tTime step: " << simulationParameters.m_timeStep << std::endl;
		stream << "\tGeometry wait time: " << simulationParameters.m_geometryWaitTime << std::endl;
		stream << "\tGeometry configurations: " << simulationParameters.m_geometryConfigurations << std::endl;
		stream << "\tEnergy output file path: " << simulationParameters.m_energyOutputFilePath << std::endl;
		stream << "\tGeometry output file path: " << simulationParameters.m_geometryOutputFilePath << std::endl;
		stream << "\tEnergy wait time: " << simulationParameters.m_energyWaitTime << std::endl;
		stream << "\tEnergy configurations: " << simulationParameters.m_energyConfigurations << std::endl;
		stream << "\tStatus wait time: " << simulationParameters.m_statusWaitTime << std::endl;
		stream << "\tEquilibrium time: " << simulationParameters.m_equilibriumTime << std::endl;
		stream << "\tEquilibrium rate: " << simulationParameters.m_equilibriumRate << std::endl;
		stream << "\tRandom seed: " << simulationParameters.m_randomSeed << std::endl;
		stream << "\tEnergy write digits: " << simulationParameters.m_energyWriteDigits << std::endl;
		stream << "\tEnergy write characters: " << simulationParameters.m_energyWriteChars << std::endl;
		stream << "\tGeometry write digits: " << simulationParameters.m_geometryWriteDigits << std::endl;
		stream << "\tGeometry write characters: " << simulationParameters.m_geometryWriteChars << std::endl;
		stream << "\tTime write digits: " << simulationParameters.m_timeWriteDigits << std::endl;
		stream << "\tTime write characters: " << simulationParameters.m_timeWriteChars << std::endl;
		stream << "]";
		return stream;
	}

	void SimulationParameters::SetDefaultParameters() {
		m_desiredTemperature = 298.15;
		m_desiredPressure = 1.0;
		m_boundarySpring = 250.0;
		m_boundary = 10.0;
		m_boundaryType = "sphere";
		m_origin = math::Vec3();
		m_totalTime = 0.5;
		m_totalConfigurations = 1000;
		m_timeStep = 0.0005;
		m_geometryWaitTime = 0.001;
		m_geometryConfigurations = 1;
		m_energyOutputFilePath = "energy.dat";
		m_geometryOutputFilePath = "geometry.xyz";
		m_energyWaitTime = 0.001;
		m_energyConfigurations = 1;
		m_statusWaitTime = 5.0;
		m_equilibriumTime = 0.0;
		m_equilibriumRate = 2.0;
		m_randomSeed = rand();
		m_energyWriteDigits = 3;
		m_energyWriteChars = 10;
		m_geometryWriteDigits = 3;
		m_geometryWriteChars = 7;
		m_timeWriteDigits = 4;
		m_timeWriteChars = 7;
	}

	void SimulationParameters::ResolveParameter(const String &key, const String &value) {
		if (key.find("desired-temperature") != String::npos) { m_desiredTemperature = utils::ToDouble(value); }
		if (key.find("desired-pressure") != String::npos) { m_desiredPressure = utils::ToDouble(value); }
		if (key.find("boundary-spring") != String::npos) { m_boundarySpring = utils::ToDouble(value); }
		if (key.find("boundary") != String::npos) { m_boundary = utils::ToDouble(value); }
		if (key.find("boundary-type") != String::npos) { m_boundaryType = value; }
		if (key.find("origin") != String::npos) {
			std::vector<String> valueTokens = utils::SplitString(value, ',');

			if (valueTokens.size() >= 3) {
				m_origin = math::Vec3(utils::ToDouble(valueTokens[0]), utils::ToDouble(valueTokens[1]), utils::ToDouble(valueTokens[2]));
			}
		}
		if (key.find("total-time") != String::npos) { m_totalTime = utils::ToDouble(value); }
		if (key.find("total-configurations") != String::npos) { m_totalConfigurations = utils::NextInt(value); }
		if (key.find("time-step") != String::npos) { m_timeStep = utils::ToDouble(value); }
		if (key.find("geometry-wait-time") != String::npos) { m_geometryWaitTime = utils::ToDouble(value); }
		if (key.find("geometry-configurations") != String::npos) { m_geometryConfigurations = utils::NextInt(value); }
		if (key.find("energy-output-file-path") != String::npos) { m_energyOutputFilePath = value; }
		if (key.find("geometry-output-file-path") != String::npos) { m_geometryOutputFilePath = value; }
		if (key.find("energy-wait-time") != String::npos) { m_energyWaitTime = utils::ToDouble(value); }
		if (key.find("energy-configurations") != String::npos) { m_energyConfigurations = utils::NextInt(value); }
		if (key.find("status-wait-time") != String::npos) { m_statusWaitTime = utils::ToDouble(value); }
		if (key.find("equilibrium-time") != String::npos) { m_equilibriumTime = utils::ToDouble(value); }
		if (key.find("equilibrium-rate") != String::npos) { m_equilibriumRate = utils::ToDouble(value); }
		if (key.find("random-seed") != String::npos) { m_randomSeed = utils::NextInt(value); }
		if (key.find("energy-write-digits") != String::npos) { m_energyWriteDigits = utils::NextInt(value); }
		if (key.find("energy-write-chars") != String::npos) { m_energyWriteChars = utils::NextInt(value); }
		if (key.find("geometry-write-digits") != String::npos) { m_geometryWriteDigits = utils::NextInt(value); }
		if (key.find("geometry-write-chars") != String::npos) { m_geometryWriteChars = utils::NextInt(value); }
		if (key.find("time-write-digits") != String::npos) { m_timeWriteDigits = utils::NextInt(value); }
		if (key.find("time-write-chars") != String::npos) { m_timeWriteChars = utils::NextInt(value); }
	}

}