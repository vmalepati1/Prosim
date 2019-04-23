#pragma once

#include "Math/PSMath.h"
#include "Utils/String.h"

namespace classical {

	double GetR2ij(const math::Vec3 &coordsi, const math::Vec3 &coordsj);
	double GetRij(const math::Vec3 &coordsi, const math::Vec3 &coordsj);
	math::Vec3 GetUij(const math::Vec3 &coordsi, const math::Vec3 &coordsj, double rij = -1);
	double GetUdp(const math::Vec3 &uveci, const math::Vec3 &uvecj);
	math::Vec3 GetUcp(const math::Vec3 &uveci, const math::Vec3 &uvecj);
	math::Vec3 GetCp(const math::Vec3 &uveci, const math::Vec3 &uvecj);
	double GetAijk(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, double rij = -1, double rjk = -1);
	double GetTijkl(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, const math::Vec3 &coordsl, double rij = -1, double rjk = -1, double rkl = -1);
	double GetOijkl(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, const math::Vec3 &coordsl, double rki = -1, double rkj = -1, double rkl = -1);
	double GetVolume(double bound, const String &boundType);

}