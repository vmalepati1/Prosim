#include "MolecularDynamics.h"

#include <chrono>
#include <random>

#include "Constants.h"
#include "FileIO.h"

namespace classical {

	MolecularDynamics::MolecularDynamics(Molecule *molecule, const SimulationParameters &parameters)
		: Simulation(molecule, parameters) {

	}

	void MolecularDynamics::Run() {
		OpenOutputFiles();
		InitializeVelocities();
		m_molecule->CalculateEnergy();
		m_molecule->CalculateGradient();
		UpdateAccelerations();
		CheckPrint(0.0, true);
		UpdateVelocities(0.5 * m_parameters.GetTimeStep());

		while (m_currentTime < m_parameters.GetTotalTime()) {
			UpdatePositions(m_parameters.GetTimeStep());
			m_molecule->CalculateGradient();
			UpdateAccelerations();
			UpdateVelocities(m_parameters.GetTimeStep());
			m_molecule->CalculateEnergy("leapfrog");

			if (m_currentTime < m_parameters.GetEquilibriumTime()) {
				EquilibrateTemperature();
			}

			CheckPrint(m_parameters.GetTimeStep());
			m_currentTime += m_parameters.GetTimeStep();
		}

		CheckPrint(m_parameters.GetTimeStep());
		CloseOutputFiles();
	}

	void MolecularDynamics::OpenOutputFiles() {
		m_energyFile.open(m_parameters.GetEnergyOutputFilePath());
		m_geometryFile.open(m_parameters.GetGeometryOutputFilePath());

		WriteEnergyHeader();

		m_lastTime = std::chrono::duration_cast<std::chrono::duration<double>>(
			std::chrono::system_clock::now().time_since_epoch()).count();

		m_eTime = 10E-10;
		m_gTime = 10E-10;
	}

	void MolecularDynamics::WriteGeometry() {
		char comment[20];
		
		snprintf(comment, 20, "%.4f ps", m_currentTime);

		m_geometryFile << GetCoordsXYZString(m_molecule->m_atoms, String(comment), m_parameters.GetGeometryWriteChars(), m_parameters.GetGeometryWriteDigits());
	}

	void MolecularDynamics::WriteEnergyTerms(int totalFloatChars, int decimalChars, char printType) {
		Molecule *m = m_molecule;

		std::vector<double> energyTerms = {
			m->m_eKinetic, m->m_ePotential, m->m_eNonBonded, m->m_eBonded, m->m_eBound,
			m->m_eVDW, m->m_eElst, m->m_eBonds, m->m_eAngles, m->m_eTorsions, m->m_eOutOfPlanes
		};

		for (double term : energyTerms) {
			WriteValue(totalFloatChars, decimalChars, term, printType);
		}
	}

	void MolecularDynamics::WriteEnergy() {
		WriteValue(m_parameters.GetTimeWriteChars(), m_parameters.GetTimeWriteDigits(), m_currentTime, 'f', 0);

		WriteValue(m_parameters.GetEnergyWriteChars() + 2, m_parameters.GetEnergyWriteDigits() + 2, m_molecule->m_eTotal, 'e');
		WriteEnergyTerms(m_parameters.GetEnergyWriteChars(), m_parameters.GetEnergyWriteDigits(), 'e');
		m_energyFile.write("\n", 1);
	}


	void MolecularDynamics::WriteEnergyHeader() {
		m_energyFile << utils::StringWithFormat("#\n# INPUTFILE %s", m_parameters.GetFilePath().c_str()).c_str();
		m_energyFile << utils::StringWithFormat("\n# ENERGYOUT %s", m_parameters.GetEnergyOutputFilePath().c_str()).c_str();
		m_energyFile << utils::StringWithFormat("\n# GEOMOUT %s", m_parameters.GetGeometryOutputFilePath().c_str()).c_str();
		m_energyFile << utils::StringWithFormat("\n# RANDOMSEED %i", m_parameters.GetRandomSeed());
		m_energyFile << utils::StringWithFormat("\n# DESIREDTEMPERATURE %.6f K", m_parameters.GetDesiredTemperature());
		m_energyFile << utils::StringWithFormat("\n# BOUNDARY %.6f A", m_molecule->m_boundary);
		m_energyFile << utils::StringWithFormat("\n# BOUNDARYSPRING %.6f kcal/(mol*A^2)", m_molecule->m_kBox);
		m_energyFile << utils::StringWithFormat("\n# BOUNDARYTYPE %s", m_molecule->m_boundaryType.c_str());
		m_energyFile << utils::StringWithFormat("\n# STATUSWAITTIME %.6f s", m_parameters.GetStatusWaitTime());
		m_energyFile << utils::StringWithFormat("\n# ENERGYWAITTIME %.6f ps", m_parameters.GetEnergyWaitTime());
		m_energyFile << utils::StringWithFormat("\n# GEOMWAITTIME %.6f ps", m_parameters.GetGeometryWaitTime());
		m_energyFile << utils::StringWithFormat("\n# TOTALTIME %.6f ps", m_parameters.GetTotalTime());
		m_energyFile << utils::StringWithFormat("\n# TIMESTEP %.6f ps", m_parameters.GetTimeStep());
		m_energyFile << utils::StringWithFormat("\n# EQTIME %.6f ps", m_parameters.GetEquilibriumTime());
		m_energyFile << utils::StringWithFormat("\n# EQRATE %.6f p", m_parameters.GetEquilibriumRate());
		m_energyFile << "\n#\n# -- ENERGY DATA --\n#";
		m_energyFile << "\n# energy terms [kcal/mol]\n#  time      e_total      ";
		m_energyFile << "e_kin      e_pot  e_nonbond   e_bonded e_boundary      ";
		m_energyFile << "e_vdw     e_elst     e_bond    e_angle     e_tors      ";
		m_energyFile << "e_oop\n";
	}

	void MolecularDynamics::PrintStatus() {
		std::cout << utils::StringWithFormat("%.*f/%.*f ps", m_parameters.GetTimeWriteDigits(), m_currentTime, m_parameters.GetTimeWriteDigits(), m_parameters.GetTotalTime());

		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer, 80, "%H:%M:%S", timeinfo);

		std::cout << "as of " << buffer << std::endl;

		FlushBuffers();
	}

	void MolecularDynamics::InitializeVelocities() {
		if (m_parameters.GetDesiredTemperature()) {
			m_eTemperature = m_parameters.GetDesiredTemperature();

			double sigmaBase = sqrt(2.0 * GAS_CONSTANT * m_parameters.GetDesiredTemperature() / 3);

			for (Atom *atom : m_molecule->m_atoms) {
				double sigma = sigmaBase * pow(atom->mass, -0.5);

				for (int j = 0; j < 3; j++) {
					std::random_device rd{};
					std::mt19937 gen{ rd() };

					std::normal_distribution<double> distribution(0.0, sigma);

					atom->velocity[j] = distribution(gen);
				}

			}

			m_molecule->CalculateEnergy();
			m_molecule->CalculateTemperature();

			double vScale = sqrt(m_parameters.GetDesiredTemperature() / m_molecule->m_temperature);

			for (Atom *atom : m_molecule->m_atoms) {
				atom->velocity *= vScale;
			}
		}
	}

	void MolecularDynamics::EquilibrateTemperature() {
		double timeScale = m_parameters.GetTimeStep() / std::max(m_parameters.GetTimeStep(), m_parameters.GetEquilibriumRate());
		double timeWeight = 10.0 * m_parameters.GetTimeStep();

		m_eTemperature = (m_eTemperature + timeWeight * m_molecule->GetMemberTemperature()) / (1.0 + timeWeight);

		double velocityScale = 1.0 + timeScale * (sqrt(m_parameters.GetDesiredTemperature() / m_eTemperature) - 1.0);

		for (int i = 0; i < m_molecule->m_nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				m_molecule->m_atoms[i]->velocity[j] *= velocityScale;
			}
		}
	}

	void MolecularDynamics::UpdateAccelerations() {
		for (int i = 0; i < m_molecule->m_nAtoms; i++) {
			double mass = m_molecule->m_atoms[i]->mass;

			for (int j = 0; j < 3; j++) {
				m_molecule->m_atoms[i]->previousAcceleration[j] = m_molecule->m_atoms[i]->acceleration[j];
				m_molecule->m_atoms[i]->acceleration[j] = (-ACCELERATION_CONVERSION * m_molecule->m_gTotal[i][j] / mass);
			}
		}
	}

	void MolecularDynamics::UpdateVelocities(double deltaTime) {
		for (int i = 0; i < m_molecule->m_nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				m_molecule->m_atoms[i]->previousVelocity[j] = m_molecule->m_atoms[i]->velocity[j];
				m_molecule->m_atoms[i]->velocity[j] += m_molecule->m_atoms[i]->acceleration[j] * deltaTime;
			}
		}
	}

	void MolecularDynamics::UpdatePositions(double deltaTime) {
		for (int i = 0; i < m_molecule->m_nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				m_molecule->m_atoms[i]->position[j] += m_molecule->m_atoms[i]->velocity[j] * deltaTime;
			}
		}

		m_molecule->UpdateInternals();
	}

	void MolecularDynamics::CheckPrint(double timeStep, bool printAll) {
		if (printAll || m_eTime >= m_parameters.GetEnergyWaitTime()) {
			WriteEnergy();
			m_eTime = 1.0E-10;
		}

		if (printAll || m_gTime >= m_parameters.GetGeometryWaitTime()) {
			WriteGeometry();
			m_gTime = 1.0E-10;
		}

		double currentTime = std::chrono::duration_cast<std::chrono::duration<double>>(
			std::chrono::system_clock::now().time_since_epoch()).count();

		if (printAll || (currentTime - m_lastTime) > m_parameters.GetStatusWaitTime()) {
			PrintStatus();
			m_lastTime = currentTime;
		}

		m_eTime += timeStep;
		m_gTime += timeStep;
	}

}