#include "FileIO.h"

namespace classical {

	String GetCoordsXYZString(const std::vector<Atom *> &atoms, const String &comment, int totalChars, int decimalChars) {
		String string = utils::StringWithFormat("%i\n%s\n", atoms.size(), comment.c_str());

		for (Atom *atom : atoms) {
			string.append(utils::StringWithFormat("%-2s", atom->element.c_str()));

			for (int j = 0; j < 3; j++) {
				string.append(utils::StringWithFormat(" %*.*f", totalChars, decimalChars, atom->position[j]));
			}

			string.append("\n");
		}

		return string;
	}

}