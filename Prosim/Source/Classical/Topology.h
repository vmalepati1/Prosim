#pragma once

#include <vector>

#include "Angle.h"
#include "Atom.h"
#include "Bond.h"
#include "OutOfPlane.h"
#include "Torsion.h"

namespace classical {

	void CalculateBondGraph(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, ForceField *forceField);
	void CalculateBondGraphFromBonds(int natoms, std::map<int, std::map<int, double>> &bondGraph, const std::vector<Bond *> &bonds);
	void CalculateBonds(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Bond *> &bonds, ForceField *forceField);
	void CalculateAngles(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Angle *> &angles, ForceField *forceField);
	void CalculateTorsions(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Torsion *> &torsions, ForceField *forceField);
	void CalculateOutOfPlanes(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<OutOfPlane *> &outOfPlanes, ForceField *forceField);
	void CalculateNonInts(const std::vector<Bond *> &bonds, const std::vector<Angle *> &angles, const std::vector<Torsion *> &torsions, std::vector<int> &nonInts);
	void UpdateBonds(std::vector<Bond *> &bonds, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph);
	void UpdateAngles(std::vector<Angle *> &angles, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph);
	void UpdateTorsions(std::vector<Torsion *> &torsions, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph);
	void UpdateOutOfPlanes(std::vector<OutOfPlane *> &outOfPlanes, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph);

}