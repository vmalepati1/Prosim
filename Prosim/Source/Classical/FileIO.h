#pragma once

#include <vector>

#include "Atom.h"

#include "Utils/String.h"

namespace classical {

	String GetCoordsXYZString(const std::vector<Atom *> &atoms, const String &comment, int totalChars = 12, int decimalChars = 6);

}