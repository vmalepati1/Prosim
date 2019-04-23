#include "Topology.h"

#include <algorithm>

#include "Constants.h"
#include "Geometry.h"

#include "Utils/IterationTools.h"

namespace classical {

	void CalculateBondGraph(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, ForceField *forceField) {
		int natoms = atoms.size();

		utils::IterationMatrix matrix(utils::CombinationsNR(natoms, 2), 2);

		utils::CombinationKN(matrix, 2, natoms);

#ifdef PS_OPTIMIZED
#pragma omp parallel for
#endif
		for (int n = 0; n < utils::CombinationsNR(natoms, 2); n++) {
			int i = matrix(n, 0);
			int j = matrix(n, 1);

			Atom *atom1 = atoms[i];
			Atom *atom2 = atoms[j];

			double covalentRadius1 = forceField->GetCovalentRadius(atom1->element);
			double covalentRadius2 = forceField->GetCovalentRadius(atom2->element);

			double threshold = BOND_THRESHOLD * (covalentRadius1 + covalentRadius2);
			double distance = GetRij(atom1->position, atom2->position);
			double distance2 = distance * distance;

			if (distance2 < threshold * threshold) {
				bondGraph[i][j] = bondGraph[j][i] = distance;
			}
		}
	}

	void CalculateBondGraphFromBonds(int natoms, std::map<int, std::map<int, double>> &bondGraph, const std::vector<Bond *> &bonds) {
		for (auto &bond : bonds) {
			bondGraph[bond->atom1][bond->atom2] = bondGraph[bond->atom2][bond->atom1] = bond->distance;
		}
	}

	void CalculateBonds(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Bond *> &bonds, ForceField *forceField) {
		bonds.clear();

		int i = 0;

		for (auto &atom1 : atoms) {
			for (auto pair : bondGraph[i]) {
				int j = pair.first;
				if (i > j) continue;

				Atom *atom2 = atoms[j];
				double distance = bondGraph[i][j];

				double equilibriumDistance = forceField->GetBondEquilibriumLength(atom1->type, atom2->type);
				double springConstant = forceField->GetBondSpringConstant(atom1->type, atom2->type);

				if (springConstant) {
					bonds.push_back(new Bond(i, j, distance, equilibriumDistance, springConstant));
				}
			}

			i++;
		}

		std::sort(std::begin(bonds), std::end(bonds), [](Bond *a, Bond *b) {return a->atom1 < b->atom1; });
	}

	void CalculateAngles(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Angle *> &angles, ForceField *forceField) {
		angles.clear();

		int j = 0;

		for (auto &atom2 : atoms) {
			std::vector<int> row;
			std::vector<std::vector<int>> indices;

			for (std::map<int, double>::iterator it = bondGraph[j].begin(); it != bondGraph[j].end(); ++it) {
				row.push_back(it->first);
			}

			if (row.size()) {
				int *data = new int[row.size()];

				utils::CombinationArray(&row[0], data, indices, 0, row.size() - 1, 0, 2);

				delete[] data;
			}

			for (auto &v : indices) {
				int i = v[0];
				int k = v[1];

				if (i > k) continue;

				Atom *atom1 = atoms[i];
				Atom *atom3 = atoms[k];

				double aijk = GetAijk(atom1->position, atom2->position, atom3->position);
				double angleEquilibriumDegrees = forceField->GetAngleEquilibriumDegrees(atom1->type, atom2->type, atom3->type);
				double angleSpringConstant = forceField->GetAngleSpringConstant(atom1->type, atom2->type, atom3->type);

				if (angleSpringConstant) {
					angles.push_back(new Angle(i, j, k, aijk, angleEquilibriumDegrees, angleSpringConstant));
				}
			}

			j++;
		}

		std::sort(std::begin(angles), std::end(angles), [](Angle *a, Angle *b) {return a->atom1 < b->atom1; });
	}

	void CalculateTorsions(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<Torsion *> &torsions, ForceField *forceField) {
		torsions.clear();

		int j = 0;

		for (auto &atom2 : atoms) {
			std::vector<int> row;

			for (std::map<int, double>::iterator it = bondGraph[j].begin(); it != bondGraph[j].end(); ++it) {
				row.push_back(it->first);
			}

			utils::IterationMatrix matrix(utils::PermutationsNR(row.size(), 2), 2);
			utils::PermutationPairsVector(row, matrix);

			for (int n = 0; n < utils::PermutationsNR(row.size(), 2); n++) {
				int i = matrix(n, 0);
				int k = matrix(n, 1);

				if (j > k) continue;

				Atom *atom1 = atoms[i];
				Atom *atom3 = atoms[k];

				for (auto pair : bondGraph[k]) {
					int l = pair.first;

					if (l == i || l == j) continue;

					Atom *atom4 = atoms[l];

					double tijkl = GetTijkl(atom1->position, atom2->position, atom3->position, atom4->position);

					std::vector<std::tuple<double, double, int, int>> torsionParameters = 
						forceField->GetTorsionParameters(atom1->type, atom2->type, atom3->type, atom4->type);
					
					for (auto &tuple : torsionParameters) {
						double vn = std::get<0>(tuple);
						double gamma = std::get<1>(tuple);
						int nfold = std::get<2>(tuple);
						int paths = std::get<3>(tuple);

						if (vn) {
							torsions.push_back(new Torsion(i, j, k, l, tijkl, vn, gamma, nfold, paths));
						}
					}
				}
			}

			j++;
		}

		std::sort(std::begin(torsions), std::end(torsions), [](Torsion *a, Torsion *b) {return a->atom1 < b->atom1; });
	}

	void CalculateOutOfPlanes(const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph, std::vector<OutOfPlane *> &outOfPlanes, ForceField *forceField) {
		outOfPlanes.clear();

		int k = 0;

		for (auto &atom : atoms) {
			std::vector<int> row;
			std::vector<std::vector<int>> indices;

			for (std::map<int, double>::iterator it = bondGraph[k].begin(); it != bondGraph[k].end(); ++it) {
				row.push_back(it->first);
			}

			if (row.size()) {
				int *data = new int[row.size()];

				utils::CombinationArray(&row[0], data, indices, 0, row.size() - 1, 0, 3);

				delete[] data;
			}

			for (auto &v : indices) {
				int i = v[0];
				int j = v[1];
				int l = v[2];

				Atom *atom1 = atoms[i];
				Atom *atom2 = atoms[j];
				Atom *atom3 = atoms[l];

				std::vector<std::tuple<int, int, int, int>> combos = {
					{ std::min(i, j), std::max(i, j), k, l },
					{ std::min(j, l), std::max(j, l), k, i },
					{ std::min(l, i), std::max(l, i), k, j },
				};

				for (auto combo : combos) {
					int x = std::get<0>(combo);
					int y = std::get<1>(combo);
					int z = std::get<2>(combo);
					int w = std::get<3>(combo);

					double oijkl = GetOijkl(atoms[x]->position, atoms[y]->position, atoms[z]->position, atoms[w]->position);
					double vn = forceField->GetOutOfPlaneHeightBarrier(atoms[x]->type, atoms[y]->type, atoms[z]->type, atoms[w]->type);

					if (vn) {
						outOfPlanes.push_back(new OutOfPlane(x, y, z, w, oijkl, vn));
					}
				}
			}

			k++;
		}

		std::sort(std::begin(outOfPlanes), std::end(outOfPlanes), [](OutOfPlane *a, OutOfPlane *b) {return a->atom1 < b->atom1; });
	}

	void CalculateNonInts(const std::vector<Bond *> &bonds, const std::vector<Angle *> &angles, const std::vector<Torsion *> &torsions, std::vector<int> &nonInts) {
		nonInts.clear();

		for (Bond *bond : bonds) {
			nonInts.push_back(bond->atom1);
			nonInts.push_back(bond->atom2);
			nonInts.push_back(bond->atom2);
			nonInts.push_back(bond->atom1);
		}

		for (Angle *angle : angles) {
			nonInts.push_back(angle->atom1);
			nonInts.push_back(angle->atom3);
			nonInts.push_back(angle->atom3);
			nonInts.push_back(angle->atom1);
		}

		for (Torsion *torsion : torsions) {
			nonInts.push_back(torsion->atom1);
			nonInts.push_back(torsion->atom4);
			nonInts.push_back(torsion->atom4);
			nonInts.push_back(torsion->atom1);
		}

		// TODO: implement sorting. std::sort(nonInts.begin(), nonInts.end());
	}

	void UpdateBonds(std::vector<Bond *> &bonds, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph) {
		for (Bond *bond : bonds) {
			math::Vec3 position1 = atoms[bond->atom1]->position;
			math::Vec3 position2 = atoms[bond->atom2]->position;

			bond->distance = GetRij(position1, position2);

			bondGraph[bond->atom1][bond->atom2] = bond->distance;
			bondGraph[bond->atom2][bond->atom1] = bond->distance;
		}
	}

	void UpdateAngles(std::vector<Angle *> &angles, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph) {
		for (Angle *angle : angles) {
			double r12 = bondGraph[angle->atom1][angle->atom2];
			double r23 = bondGraph[angle->atom2][angle->atom3];
			math::Vec3 position1 = atoms[angle->atom1]->position;
			math::Vec3 position2 = atoms[angle->atom2]->position;
			math::Vec3 position3 = atoms[angle->atom3]->position;

			angle->degrees = GetAijk(position1, position2, position3, r12, r23);
		}
	}

	void UpdateTorsions(std::vector<Torsion *> &torsions, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph) {
		for (Torsion *torsion : torsions) {
			double r12 = bondGraph[torsion->atom1][torsion->atom2];
			double r23 = bondGraph[torsion->atom2][torsion->atom3];
			double r34 = bondGraph[torsion->atom3][torsion->atom4];
			math::Vec3 position1 = atoms[torsion->atom1]->position;
			math::Vec3 position2 = atoms[torsion->atom2]->position;
			math::Vec3 position3 = atoms[torsion->atom3]->position;
			math::Vec3 position4 = atoms[torsion->atom4]->position;

			torsion->degrees = GetTijkl(position1, position2, position3, position4, r12, r23, r34);
		}
	}

	void UpdateOutOfPlanes(std::vector<OutOfPlane *> &outOfPlanes, const std::vector<Atom *> &atoms, std::map<int, std::map<int, double>> &bondGraph) {
		for (OutOfPlane *outOfPlane : outOfPlanes) {
			double r31 = bondGraph[outOfPlane->atom3][outOfPlane->atom1];
			double r32 = bondGraph[outOfPlane->atom3][outOfPlane->atom2];
			double r34 = bondGraph[outOfPlane->atom3][outOfPlane->atom4];
			math::Vec3 position1 = atoms[outOfPlane->atom1]->position;
			math::Vec3 position2 = atoms[outOfPlane->atom2]->position;
			math::Vec3 position3 = atoms[outOfPlane->atom3]->position;
			math::Vec3 position4 = atoms[outOfPlane->atom4]->position;

			outOfPlane->degrees = GetOijkl(position1, position2, position3, position4, r31, r32, r34);
		}
	}

}