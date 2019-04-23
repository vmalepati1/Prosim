#include <chrono>
#include <CL/cl.h>
#include <iostream>

#include "Source/Classical/SimulationParameters.h"
#include "Source/Classical/ForceField.h"
#include "Source/Classical/Atom.h"
#include "Source/Classical/Molecule.h"
#include "Source/Classical/PQRMolecule.h"
#include "Source/Classical/MolecularDynamics.h"
#include "Source/Classical/Utils/IterationTools.h"

using namespace classical;
using namespace classical::math;

int main() {
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	SimulationParameters simulationParameters;

	ForceField forceField;

	PQRMolecule molecule("Tests/Ethane.pqr", &forceField, true);

	MolecularDynamics md(&molecule, simulationParameters);

	md.Run();

	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	std::cout << "Time taken for simulation (microseconds): " << duration << std::endl;

	system("PAUSE");

	return 0;
}