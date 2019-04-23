#pragma once

#include "Utils/String.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

namespace classical {

	class ForceField {
	public:
		ForceField();
		ForceField(const String &filePath);

		inline double GetAtomicMass(const String &element) { 
			if (m_atomicMasses.find(element) != m_atomicMasses.end()) return m_atomicMasses.find(element)->second;

			std::cout << "Could not find atomic mass for element " << element << std::endl;
			return -1;
		}

		inline double GetCovalentRadius(const String &element) { 
			if (m_covalentRadii.find(element) != m_covalentRadii.end()) return m_covalentRadii.find(element)->second;

			std::cout << "Could not find covalent radius for element " << element << std::endl;
			return -1;
		}

		inline double GetVanDerWaalsRadius(const String &type) { 
			if (m_vanDerWaalsParameters.find(type) != m_vanDerWaalsParameters.end()) return m_vanDerWaalsParameters.find(type)->second.first;

			std::cout << "Could not find Van Der Waals radius for atom type " << type << std::endl;
			return -1;
		}

		inline double GetVanDerWaalsAttractionMagnitude(const String &type) { 
			if (m_vanDerWaalsParameters.find(type) != m_vanDerWaalsParameters.end()) return m_vanDerWaalsParameters.find(type)->second.second;

			std::cout << "Could not find Van Der Waals attraction magnitude for atom type " << type << std::endl;
			return -1;
		}

		inline double GetBondSpringConstant(const String &type1, const String &type2) { 
			if (m_bondLengthParameters.find(std::make_pair(type1, type2)) != m_bondLengthParameters.end()) {
				return m_bondLengthParameters.find(std::make_pair(type1, type2))->second.first;
			}
			
			if (m_bondLengthParameters.find(std::make_pair(type2, type1)) != m_bondLengthParameters.end()) {
				return m_bondLengthParameters.find(std::make_pair(type2, type1))->second.first;
			}

			std::cout << "Could not find bond spring constant for atom type pair (" << type1 << ", " << type2 << ")" << std::endl;
			return -1;
		}

		inline double GetBondEquilibriumLength(const String &type1, const String &type2) { 
			if (m_bondLengthParameters.find(std::make_pair(type1, type2)) != m_bondLengthParameters.end()) {
				return m_bondLengthParameters.find(std::make_pair(type1, type2))->second.second;
			}
			
			if (m_bondLengthParameters.find(std::make_pair(type2, type1)) != m_bondLengthParameters.end()) {
				return m_bondLengthParameters.find(std::make_pair(type2, type1))->second.second;
			}

			std::cout << "Could not find bond equilibrium length for atom type pair (" << type1 << ", " << type2  << ")" << std::endl;
			return -1;
		}

		inline double GetAngleSpringConstant(const String &type1, const String &type2, const String &type3) {
			if (m_bondAngleParameters.find(std::make_tuple(type1, type2, type3)) != m_bondAngleParameters.end()) {
				return m_bondAngleParameters.find(std::make_tuple(type1, type2, type3))->second.first;
			}

			if (m_bondAngleParameters.find(std::make_tuple(type3, type2, type1)) != m_bondAngleParameters.end()) {
				return m_bondAngleParameters.find(std::make_tuple(type3, type2, type1))->second.first;
			}

			std::cout << "Could not find angle spring constant for atom type triplet (" << type1 << ", " << type2 << ", " << type3 << ")" << std::endl;
			return -1;
		}

		inline double GetAngleEquilibriumDegrees(const String &type1, const String &type2, const String &type3) {
			if (m_bondAngleParameters.find(std::make_tuple(type1, type2, type3)) != m_bondAngleParameters.end()) {
				return m_bondAngleParameters.find(std::make_tuple(type1, type2, type3))->second.second;
			}

			if (m_bondAngleParameters.find(std::make_tuple(type3, type2, type1)) != m_bondAngleParameters.end()) {
				return m_bondAngleParameters.find(std::make_tuple(type3, type2, type1))->second.second;
			}

			std::cout << "Could not find angle equilibrium degrees for atom type triplet (" << type1 << ", " << type2 << ", " << type3 << ")" << std::endl;
			return -1;
		}

		std::vector<std::tuple<double, double, int, int>> GetTorsionParameters(const String &type1, const String &type2, const String &type3, const String &type4) {
			std::vector<std::tuple<double, double, int, int>> torsionParameters;

			if (m_torsion23Parameters.find(std::make_pair(type2, type3)) != m_torsion23Parameters.end()) {
				torsionParameters.push_back(m_torsion23Parameters.find(std::make_pair(type2, type3))->second);
			} else if (m_torsion23Parameters.find(std::make_pair(type3, type2)) != m_torsion23Parameters.end()) {
				torsionParameters.push_back(m_torsion23Parameters.find(std::make_pair(type3, type2))->second);
			}
			else {
				std::cout << "Could not find torsion parameters for atom central atom pair (" << type2 << ", " << type3 << ")" << std::endl;
			}

			std::map<std::tuple<String, String, String, String>, std::vector<std::tuple<double, double, int, int>>>::iterator foundTorsion1234Iterator = 
				m_torsion1234Parameters.find(std::make_tuple(type1, type2, type3, type4));

			if (foundTorsion1234Iterator != m_torsion1234Parameters.end()) {
				torsionParameters.insert(torsionParameters.end(), foundTorsion1234Iterator->second.begin(), foundTorsion1234Iterator->second.end());
			}

			foundTorsion1234Iterator = m_torsion1234Parameters.find(std::make_tuple(type4, type3, type2, type1));

			if (foundTorsion1234Iterator != m_torsion1234Parameters.end()) {
				torsionParameters.insert(torsionParameters.end(), foundTorsion1234Iterator->second.begin(), foundTorsion1234Iterator->second.end());
			}

			return torsionParameters;
		}

		inline double GetOutOfPlaneHeightBarrier(const String &type1, const String &type2, const String &type3, const String &type4) {
			if (m_outOfPlane1234Parameters.find(std::make_tuple(type1, type2, type3, type4)) != m_outOfPlane1234Parameters.end()) {
				return m_outOfPlane1234Parameters.find(std::make_tuple(type1, type2, type3, type4))->second;
			}

			if (m_outOfPlane234Parameters.find(std::make_tuple(type2, type3, type4)) != m_outOfPlane234Parameters.end()) {
				return m_outOfPlane234Parameters.find(std::make_tuple(type2, type3, type4))->second;
			}

			if (m_outOfPlane34Parameters.find(std::make_pair(type3, type4)) != m_outOfPlane34Parameters.end()) {
				return m_outOfPlane34Parameters.find(std::make_pair(type3, type4))->second;
			}

			return 0.0;
		}

		inline const std::map<String, double> &GetAtomicMasses() const { return m_atomicMasses; }
		inline const std::map<String, double> &GetElementsRadii() const { return m_covalentRadii; }
		inline const std::map<String, std::pair<double, double>> &GetVanDerWaalsParameters() const { return m_vanDerWaalsParameters; }
		inline const std::map<std::pair<String, String>, std::pair<double, double>> &GetBondLengthParameters() const { return m_bondLengthParameters; }
		inline const std::map<std::tuple<String, String, String>, std::pair<double, double>> &GetBondAngleParameters() const { return m_bondAngleParameters; }
		inline const std::map<std::pair<String, String>, std::tuple<double, double, int, int>> &GetTorsion23Parameters() const { return m_torsion23Parameters; }
		inline const std::map<std::tuple<String, String, String, String>, std::vector<std::tuple<double, double, int, int>>> &GetTorsion1234Parameters() const { return m_torsion1234Parameters; }
		inline const std::map<std::pair<String, String>, double> &GetOutOfPlane34Parameters() const { return m_outOfPlane34Parameters; }
		inline const std::map<std::tuple<String, String, String>, double> &GetOutOfPlane234Parameters() const { return m_outOfPlane234Parameters; }
		inline const std::map<std::tuple<String, String, String, String>, double> &GetOutOfPlane1234Parameters() const { return m_outOfPlane1234Parameters; }

		friend std::ostream& operator<<(std::ostream &stream, const ForceField &forceField);
	private:
		void SetDefaultParameters();
		void ResolveTokens(const std::vector<String> &tokens);
	private:
		/* String: atomic symbol, Double: atomic mass. */
		std::map<String, double> m_atomicMasses;
		/* String: atomic symbol, Double: radius. */
		std::map<String, double> m_covalentRadii;
		/* String: atom type, Double 1: Van Der Waals radius ro/2 [Angstrom], Double 2: Van Der Waals attraction magnitude eps [kcal/mol]. */
		std::map<String, std::pair<double, double>> m_vanDerWaalsParameters;
		/* String 1: atom type, String 2: atom type, Double 1: bond spring constant k_b [kcal/(mol*A^2)], Double 2: bond equilibrium length r_eq [Angstrom]. */
		std::map<std::pair<String, String>, std::pair<double, double>> m_bondLengthParameters;
		/* String 1: atom type, String 2: atom type, String 3: atom type, Double 1: angle spring constant k_a [kcal/(mol*rad^2)], Double 2: bond equilibrium angle a_eq [degrees]. */
		std::map<std::tuple<String, String, String>, std::pair<double, double>> m_bondAngleParameters;
		/* String 1: atom type, String 2: atom type, Double 1: rotation barrier height vn/2 [kcal/mol], Double 2: barrier minimum offset angle gamma [degrees], Int 1: frequency of barrier n, Int 2: number of unique torsion paths paths. */
		std::map<std::pair<String, String>, std::tuple<double, double, int, int>> m_torsion23Parameters;
		/* String 1: atom type, String 2: atom type, String 3: atom type, String 4: atom type, Doubles/Ints: same as torsion23 parameters (variable length). */
		std::map<std::tuple<String, String, String, String>, std::vector<std::tuple<double, double, int, int>>> m_torsion1234Parameters;
		/* String 1: atom type, String 2: atom type, Double: rotation barrier height vn/2 [kcal/mol]. */
		std::map<std::pair<String, String>, double> m_outOfPlane34Parameters;
		/* String 1: atom type, String 2: atom type, String 3: atom type, Double: rotation barrier height vn/2 [kcal/mol]. */
		std::map<std::tuple<String, String, String>, double> m_outOfPlane234Parameters;
		/* String 1: atom type, String 2: atom type, String 3: atom type, String 4: atom type, Double: rotation barrier height vn/2 [kcal/mol]. */
		std::map<std::tuple<String, String, String, String>, double> m_outOfPlane1234Parameters;
	};

}