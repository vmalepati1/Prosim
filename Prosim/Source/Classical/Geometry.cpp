#include "Geometry.h"

#include <algorithm>
#include <math.h>

#include "Constants.h"

namespace classical {

	double GetR2ij(const math::Vec3 &coordsi, const math::Vec3 &coordsj) {
		return pow(coordsi.Distance(coordsj), 2);
	}

	double GetRij(const math::Vec3 &coordsi, const math::Vec3 &coordsj) {
		return coordsi.Distance(coordsj);
	}

	math::Vec3 GetUij(const math::Vec3 &coordsi, const math::Vec3 &coordsj, double rij) {
		
		if (rij < 0) rij = GetRij(coordsi, coordsj);

		if (!rij) return math::Vec3();

		return math::Vec3(coordsj.x - coordsi.x, coordsj.y - coordsi.y, coordsj.z - coordsi.z) / rij;
	}

	double GetUdp(const math::Vec3 &uveci, const math::Vec3 &uvecj) {
		double udp;

		udp = uveci.x * uvecj.x;
		udp += uveci.y * uvecj.y;
		udp += uveci.z * uvecj.z;

		return std::max(-1.0, std::min(1.0, udp));
	}

	math::Vec3 GetUcp(const math::Vec3 &uveci, const math::Vec3 &uvecj) {
		math::Vec3 ucp;

		double cosijk = GetUdp(uveci, uvecj);
		double sinijk = sqrt(1.0 - cosijk * cosijk);

		if (sinijk) {
			ucp.x = (uveci.y * uvecj.z - uveci.z*uvecj.y) / sinijk;
			ucp.y = (uveci.z * uvecj.x - uveci.x*uvecj.z) / sinijk;
			ucp.z = (uveci.x * uvecj.y - uveci.y*uvecj.x) / sinijk;
		}

		return ucp;
	}

	math::Vec3 GetCp(const math::Vec3 &uveci, const math::Vec3 &uvecj) {
		return uveci.Cross(uvecj);
	}

	double GetAijk(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, double rij, double rjk) {
		math::Vec3 uji = GetUij(coordsj, coordsi, rij);
		math::Vec3 ujk = GetUij(coordsj, coordsk, rjk);
		double dpjijk = GetUdp(uji, ujk);

		return RADIANS_TO_DEGREES * acos(dpjijk);
	}

	double GetTijkl(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, const math::Vec3 &coordsl, double rij, double rjk, double rkl) {
		math::Vec3 uji = GetUij(coordsj, coordsi, rij);
		math::Vec3 ujk = GetUij(coordsj, coordsk, rjk);
		math::Vec3 ukl = GetUij(coordsk, coordsl, rkl);
		math::Vec3 ujijk = GetUcp(uji, ujk);
		math::Vec3 ukjkl = GetUcp(ujk, ukl) * -1;
		double dpjijkkjkl = GetUdp(ujijk, ukjkl);
		double sign = (GetUdp(ujijk, ukl) <= 0.0 ? 1.0 : -1.0);
		return RADIANS_TO_DEGREES * sign * acos(dpjijkkjkl);
	}

	double GetOijkl(const math::Vec3 &coordsi, const math::Vec3 &coordsj, const math::Vec3 &coordsk, const math::Vec3 &coordsl, double rki, double rkj, double rkl) {
		math::Vec3 uki = GetUij(coordsk, coordsi, rki);
		math::Vec3 ukj = GetUij(coordsk, coordsj, rkj);
		math::Vec3 ukl = GetUij(coordsk, coordsl, rkl);
		math::Vec3 ukikj = GetUcp(uki, ukj);
		double dpkikjkl = GetUdp(ukikj, ukl);
		return RADIANS_TO_DEGREES * asin(dpkikjkl);
	}

	double GetVolume(double bound, const String &boundType) {
		if (boundType == "cube") {
			return 8.0 * pow(bound, 3);
		}
		else if (boundType == "sphere") {
			return 4.0 / 3.0 * M_PI * pow(bound, 3);
		}
		else {
			std::cout << "Unknown boundary type: " << boundType << std::endl;
			std::cout << "Use 'cube' or 'sphere'" << std::endl;
			return -1;
		}
	}

}