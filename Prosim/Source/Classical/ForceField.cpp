#include "ForceField.h"

#include <algorithm>
#include <fstream>
#include <sstream>

namespace classical {

	const static std::map<String, double> s_defaultAtomicMasses = {
		{ "H" , 1.00794 },
		{ "P" , 30.9738 },
		{ "He", 4.00260 },
		{ "B" , 10.8110 },
		{ "K" , 39.0983 },
		{ "Cr", 51.9961 },
		{ "Cu", 63.5460 },
		{ "Se", 78.9600 },
		{ "C" , 12.0107 },
		{ "S" , 32.0650 },
		{ "Ne", 20.1797 },
		{ "Na", 22.9898 },
		{ "Ca", 40.0780 },
		{ "Mn", 54.9380 },
		{ "Zn", 65.4090 },
		{ "Kr", 83.7980 },
		{ "O" , 15.9994 },
		{ "Cl", 35.4530 },
		{ "Ar", 39.9480 },
		{ "Mg", 24.3050 },
		{ "Sc", 44.9559 },
		{ "Fe", 55.8450 },
		{ "Ga", 69.7230 },
		{ "N" , 14.0067 },
		{ "Br", 79.9040 },
		{ "Li", 6.94100 },
		{ "Al", 26.9815 },
		{ "Ti", 47.8670 },
		{ "Co", 58.9332 },
		{ "Ge", 72.6400 },
		{ "F" , 18.9984 },
		{ "I" , 126.904 },
		{ "Be", 9.01218 },
		{ "Si", 28.0855 },
		{ "V" , 50.9415 },
		{ "Ni", 58.6934 },
		{ "As", 74.9216 },
		{ "X" , 0.00000 }
	};

	const static std::map<String, double> s_defaultCovalentRadii = {
		{ "H",  0.37 },
		{ "C",  0.77 },
		{ "O",  0.73 },
		{ "N",  0.75 },
		{ "F",  0.71 },
		{ "P",  1.10 },
		{ "S",  1.03 },
		{ "Cl", 0.99 },
		{ "Br", 1.14 },
		{ "I",  1.33 },
		{ "He", 0.30 },
		{ "Ne", 0.84 },
		{ "Ar", 1.00 },
		{ "Li", 1.02 },
		{ "Be", 0.27 },
		{ "B",  0.88 },
		{ "Na", 1.02 },
		{ "Mg", 0.72 },
		{ "Al", 1.30 },
		{ "Si", 1.18 },
		{ "K",  1.38 },
		{ "Ca", 1.00 },
		{ "Sc", 0.75 },
		{ "Ti", 0.86 },
		{ "V",  0.79 },
		{ "Cr", 0.73 },
		{ "Mn", 0.67 },
		{ "Fe", 0.61 },
		{ "Co", 0.64 },
		{ "Ni", 0.55 },
		{ "Cu", 0.46 },
		{ "Zn", 0.60 },
		{ "Ga", 1.22 },
		{ "Ge", 1.22 },
		{ "As", 1.22 },
		{ "Se", 1.17 },
		{ "Kr", 1.03 },
		{ "X",  0.00 }
	};

	const static std::map<String, std::pair<double, double>> s_defaultVanDerWaalsParameters = {
		{ "C" ,{ 1.9080, 0.0860 } },
		{ "CC",{ 1.9080, 0.0860 } },
		{ "CR",{ 1.9080, 0.0860 } },
		{ "CN",{ 1.9080, 0.0860 } },
		{ "CT",{ 1.9080, 0.1094 } },
		{ "H" ,{ 0.6000, 0.0157 } },
		{ "H3",{ 1.1870, 0.0157 } },
		{ "HA",{ 1.4590, 0.0150 } },
		{ "HP",{ 1.1000, 0.0157 } },
		{ "Li",{ 1.1370, 0.0183 } },
		{ "N" ,{ 1.8240, 0.1700 } },
		{ "O2",{ 1.6612, 0.2100 } },
		{ "OW",{ 1.7683, 0.1520 } },
		{ "SH",{ 2.0000, 0.2500 } },
		{ "NB",{ 1.8240, 0.1700 } },
		{ "N2",{ 1.8240, 0.1700 } },
		{ "Ar",{ 1.8436, 0.4466 } },
		{ "CA",{ 1.9080, 0.0860 } },
		{ "CV",{ 1.9080, 0.0860 } },
		{ "CB",{ 1.9080, 0.0860 } },
		{ "CK",{ 1.9080, 0.0860 } },
		{ "CS",{ 3.3950, 0.0000806 } },
		{ "H1",{ 1.3870, 0.0157 } },
		{ "H4",{ 1.4090, 0.0150 } },
		{ "HC",{ 1.4870, 0.0157 } },
		{ "HS",{ 0.6000, 0.0157 } },
		{ "IP",{ 1.8680, 0.00277 } },
		{ "K" ,{ 2.6580, 0.000328 } },
		{ "OH",{ 1.7210, 0.2104 } },
		{ "P" ,{ 2.1000, 0.2000 } },
		{ "N" ,{ 1.8240, 0.1700 } },
		{ "NC",{ 1.8240, 0.1700 } },
		{ "N3",{ 1.8750, 0.1700 } },
		{ "HH",{ 1.5000, 0.0000 } },
		{ "CM",{ 1.9080, 0.0860 } },
		{ "CW",{ 1.9080, 0.0860 } },
		{ "C*",{ 1.9080, 0.0860 } },
		{ "CQ",{ 1.9080, 0.0860 } },
		{ "F" ,{ 1.7500, 0.0610 } },
		{ "H2",{ 1.2870, 0.0157 } },
		{ "H5",{ 1.3590, 0.0150 } },
		{ "HO",{ 0.0001, 0.0000 } },
		{ "HW",{ 0.0001, 0.0000 } },
		{ "Li",{ 1.1370, 0.0183 } },
		{ "O" ,{ 1.6612, 0.2100 } },
		{ "OS",{ 1.6837, 0.1700 } },
		{ "S" ,{ 2.0000, 0.2500 } },
		{ "NA",{ 1.8240, 0.1700 } },
		{ "N*",{ 1.8240, 0.1700 } },
		{ "He",{ 1.5800, 0.0112 } },
		{ "X" ,{ 0.0000, 0.0000 } }
	};

	const static std::map<std::pair<String, String>, std::pair<double, double>> s_defaultBondLengthParameters = {
		{ { "C" , "CA" },{ 469.0, 1.409 } },
		{ { "C" , "CM" },{ 410.0, 1.444 } },
		{ { "C" ,  "N" },{ 490.0, 1.335 } },
		{ { "C" , "NA" },{ 418.0, 1.388 } },
		{ { "C" , "O" },{ 570.0, 1.229 } },
		{ { "C" , "OH" },{ 450.0, 1.364 } },
		{ { "C*", "CT" },{ 317.0, 1.495 } },
		{ { "C*", "HC" },{ 367.0, 1.080 } },
		{ { "CA", "CB" },{ 469.0, 1.404 } },
		{ { "CA", "CN" },{ 469.0, 1.400 } },
		{ { "CA", "H4" },{ 367.0, 1.080 } },
		{ { "CA", "N2" },{ 481.0, 1.340 } },
		{ { "CA", "NC" },{ 483.0, 1.339 } },
		{ { "CB", "CN" },{ 447.0, 1.419 } },
		{ { "CB", "NB" },{ 414.0, 1.391 } },
		{ { "CC", "CT" },{ 317.0, 1.504 } },
		{ { "CC", "CW" },{ 518.0, 1.371 } },
		{ { "CC", "NB" },{ 410.0, 1.394 } },
		{ { "CK", "N*" },{ 440.0, 1.371 } },
		{ { "CM", "CM" },{ 549.0, 1.350 } },
		{ { "CM", "H4" },{ 367.0, 1.080 } },
		{ { "CM", "HA" },{ 367.0, 1.080 } },
		{ { "CN", "NA" },{ 428.0, 1.380 } },
		{ { "CQ", "NC" },{ 502.0, 1.324 } },
		{ { "CR", "NA" },{ 477.0, 1.343 } },
		{ { "CT", "CT" },{ 310.0, 1.526 } },
		{ { "CT", "H1" },{ 340.0, 1.090 } },
		{ { "CT", "H3" },{ 340.0, 1.090 } },
		{ { "CT", "HP" },{ 340.0, 1.090 } },
		{ { "CT", "N*" },{ 337.0, 1.475 } },
		{ { "CT", "N3" },{ 367.0, 1.471 } },
		{ { "CT", "OS" },{ 320.0, 1.410 } },
		{ { "CT", "SH" },{ 237.0, 1.810 } },
		{ { "CV", "NB" },{ 410.0, 1.394 } },
		{ { "CW", "NA" },{ 410.0, 1.394 } },
		{ { "H" , "N*" },{ 434.0, 1.010 } },
		{ { "H" , "N3" },{ 434.0, 1.010 } },
		{ { "HO", "OH" },{ 553.0, 0.960 } },
		{ { "HS", "SH" },{ 274.0, 1.336 } },
		{ { "OH", "P" },{ 230.0, 1.610 } },
		{ { "OW", "HW" },{ 553.0,0.9572 } },
		{ { "HH", "HH" },{ 100.0, 0.740 } },
		{ { "C" , "CB" },{ 447.0, 1.419 } },
		{ { "C" , "CT" },{ 317.0, 1.522 } },
		{ { "C" , "N*" },{ 424.0, 1.383 } },
		{ { "C" , "NC" },{ 457.0, 1.358 } },
		{ { "C" , "O2" },{ 656.0, 1.250 } },
		{ { "C*", "CB" },{ 388.0, 1.459 } },
		{ { "C*", "CW" },{ 546.0, 1.352 } },
		{ { "CA", "CA" },{ 469.0, 1.400 } },
		{ { "CA", "CM" },{ 427.0, 1.433 } },
		{ { "CA", "CT" },{ 317.0, 1.510 } },
		{ { "CA", "HA" },{ 367.0, 1.080 } },
		{ { "CA", "NA" },{ 427.0, 1.381 } },
		{ { "CB", "CB" },{ 520.0, 1.370 } },
		{ { "CB", "N*" },{ 436.0, 1.374 } },
		{ { "CB", "NC" },{ 461.0, 1.391 } },
		{ { "CC", "CV" },{ 512.0, 1.375 } },
		{ { "CC", "NA" },{ 422.0, 1.385 } },
		{ { "CK", "H5" },{ 367.0, 1.080 } },
		{ { "CK", "NB" },{ 529.0, 1.304 } },
		{ { "CM", "CT" },{ 317.0, 1.510 } },
		{ { "CM", "H5" },{ 367.0, 1.080 } },
		{ { "CM", "N*" },{ 448.0, 1.365 } },
		{ { "CQ", "H5" },{ 367.0, 1.080 } },
		{ { "CR", "H5" },{ 367.0, 1.080 } },
		{ { "CR", "NB" },{ 488.0, 1.335 } },
		{ { "CT", "F" },{ 367.0, 1.380 } },
		{ { "CT", "H2" },{ 340.0, 1.090 } },
		{ { "CT", "HC" },{ 340.0, 1.090 } },
		{ { "CT", "N" },{ 337.0, 1.449 } },
		{ { "CT", "N2" },{ 337.0, 1.350 } },
		{ { "CT", "OH" },{ 320.0, 1.410 } },
		{ { "CT", "S" },{ 227.0, 1.810 } },
		{ { "CV", "H4" },{ 367.0, 1.080 } },
		{ { "CW", "H4" },{ 367.0, 1.080 } },
		{ { "H" , "N" },{ 434.0, 1.010 } },
		{ { "H" , "N2" },{ 434.0, 1.010 } },
		{ { "H" , "NA" },{ 434.0, 1.010 } },
		{ { "HO", "OS" },{ 553.0, 0.960 } },
		{ { "O2", "P" },{ 525.0, 1.480 } },
		{ { "OS", "P" },{ 230.0, 1.610 } },
		{ { "S" , "S" },{ 166.0, 2.038 } },
		{ { "X" , "X" },{ 0.0, 0.000 } }
	};

	const static std::map<std::tuple<String, String, String>, std::pair<double, double>> s_defaultBondAngleParameters = {
		{ {"C" , "CA", "CA"}, { 63.0, 120.00} }, 
		{ {"C" , "CB", "CB"}, { 63.0, 119.20} }, 
		{ {"C" , "CM", "CM"}, { 63.0, 120.70} }, 
		{ {"C" , "CM", "H4"}, { 35.0, 119.70} }, 
		{ {"C" , "CT", "CT"}, { 63.0, 111.10} }, 
		{ {"C" , "CT", "HC"}, { 50.0, 109.50} }, 
		{ {"C" , "CT", "N" }, { 63.0, 110.10} }, 
		{ {"C" , "N" , "CT"}, { 50.0, 121.90} }, 
		{ {"C" , "N*", "CM"}, { 70.0, 121.60} }, 
		{ {"C" , "N*", "H" }, { 30.0, 119.20} }, 
		{ {"C" , "NA", "CA"}, { 70.0, 125.20} }, 
		{ {"C" , "NC", "CA"}, { 70.0, 120.50} }, 
		{ {"C*", "CB", "CA"}, { 63.0, 134.90} }, 
		{ {"C*", "CT", "CT"}, { 63.0, 115.60} }, 
		{ {"C*", "CW", "H4"}, { 35.0, 120.00} }, 
		{ {"CA", "C" , "CA"}, { 63.0, 120.00} }, 
		{ {"CA", "CA", "CA"}, { 63.0, 120.00} }, 
		{ {"CA", "CA", "CN"}, { 63.0, 120.00} }, 
		{ {"CA", "CA", "H4"}, { 35.0, 120.00} }, 
		{ {"CA", "CB", "CB"}, { 63.0, 117.30} }, 
		{ {"CA", "CB", "NB"}, { 70.0, 132.40} }, 
		{ {"CA", "CM", "H4"}, { 35.0, 123.30} }, 
		{ {"CA", "CN", "CB"}, { 63.0, 122.70} }, 
		{ {"CA", "CT", "CT"}, { 63.0, 114.00} }, 
		{ {"CA", "N2", "CT"}, { 50.0, 123.20} }, 
		{ {"CA", "NA", "H" }, { 30.0, 118.00} }, 
		{ {"CA", "NC", "CQ"}, { 70.0, 118.60} }, 
		{ {"CB", "C" , "O" }, { 80.0, 128.80} }, 
		{ {"CB", "C*", "CW"}, { 63.0, 106.40} }, 
		{ {"CB", "CA", "HA"}, { 35.0, 120.00} }, 
		{ {"CB", "CA", "NC"}, { 70.0, 117.30} }, 
		{ {"CB", "CB", "NB"}, { 70.0, 110.40} }, 
		{ {"CB", "CN", "NA"}, { 70.0, 104.40} }, 
		{ {"CB", "N*", "CT"}, { 70.0, 125.80} }, 
		{ {"CB", "NB", "CK"}, { 70.0, 103.80} }, 
		{ {"CC", "CT", "CT"}, { 63.0, 113.10} }, 
		{ {"CC", "CV", "H4"}, { 35.0, 120.00} }, 
		{ {"CC", "CW", "H4"}, { 35.0, 120.00} }, 
		{ {"CC", "NA", "CR"}, { 70.0, 120.00} }, 
		{ {"CC", "NB", "CR"}, { 70.0, 117.00} }, 
		{ {"CK", "N*", "H" }, { 30.0, 128.80} }, 
		{ {"CM", "C" , "O" }, { 80.0, 125.30} }, 
		{ {"CM", "CA", "NC"}, { 70.0, 121.50} }, 
		{ {"CM", "CM", "H4"}, { 35.0, 119.70} }, 
		{ {"CM", "CM", "N*"}, { 70.0, 121.20} }, 
		{ {"CM", "N*", "CT"}, { 70.0, 121.20} }, 
		{ {"CN", "CA", "HA"}, { 35.0, 120.00} }, 
		{ {"CN", "NA", "H" }, { 30.0, 123.10} }, 
		{ {"CR", "NA", "H" }, { 30.0, 120.00} }, 
		{ {"CT", "C" , "N" }, { 70.0, 116.60} }, 
		{ {"CT", "C" , "O2"}, { 70.0, 117.00} }, 
		{ {"CT", "CC", "CV"}, { 70.0, 120.00} }, 
		{ {"CT", "CC", "NA"}, { 70.0, 120.00} }, 
		{ {"CT", "CT", "CT"}, { 40.0, 109.50} }, 
		{ {"CT", "CT", "H2"}, { 50.0, 109.50} }, 
		{ {"CT", "CT", "HP"}, { 50.0, 109.50} }, 
		{ {"CT", "CT", "N*"}, { 50.0, 109.50} }, 
		{ {"CT", "CT", "N3"}, { 80.0, 111.20} }, 
		{ {"CT", "CT", "OS"}, { 50.0, 109.50} }, 
		{ {"CT", "CT", "SH"}, { 50.0, 108.60} }, 
		{ {"CT", "N" , "H" }, { 30.0, 118.04} }, 
		{ {"CT", "N3", "H" }, { 50.0, 109.50} }, 
		{ {"CT", "OS", "CT"}, { 60.0, 109.50} }, 
		{ {"CT", "S" , "CT"}, { 62.0,  98.90} }, 
		{ {"CT", "SH", "SH"}, { 43.0,  96.00} }, 
		{ {"CW", "CC", "NA"}, { 70.0, 120.00} }, 
		{ {"CW", "NA", "H" }, { 30.0, 120.00} }, 
		{ {"F" , "CT", "H1"}, { 35.0, 109.50} }, 
		{ {"H" , "N2", "H" }, { 35.0, 120.00} }, 
		{ {"H1", "CT", "H1"}, { 35.0, 109.50} }, 
		{ {"H1", "CT", "N*"}, { 50.0, 109.50} }, 
		{ {"H1", "CT", "OH"}, { 50.0, 109.50} }, 
		{ {"H1", "CT", "S" }, { 50.0, 109.50} }, 
		{ {"H2", "CT", "H2"}, { 35.0, 109.50} }, 
		{ {"H2", "CT", "OS"}, { 50.0, 109.50} }, 
		{ {"H4", "CV", "NB"}, { 35.0, 120.00} }, 
		{ {"H5", "CK", "N*"}, { 35.0, 123.05} }, 
		{ {"H5", "CQ", "NC"}, { 35.0, 115.45} }, 
		{ {"H5", "CR", "NB"}, { 35.0, 120.00} }, 
		{ {"HO", "OH", "P" }, { 45.0, 108.50} }, 
		{ {"HP", "CT", "N3"}, { 50.0, 109.50} }, 
		{ {"HW", "OW", "HW"}, {100.0, 104.52} }, 
		{ {"N*", "C" , "NA"}, { 70.0, 115.40} }, 
		{ {"N*", "C" , "O" }, { 80.0, 120.90} }, 
		{ {"N*", "CK", "NB"}, { 70.0, 113.90} }, 
		{ {"N2", "CA", "N2"}, { 70.0, 120.00} }, 
		{ {"N2", "CA", "NC"}, { 70.0, 119.30} }, 
		{ {"NA", "CA", "NC"}, { 70.0, 123.30} }, 
		{ {"NA", "CR", "NB"}, { 70.0, 120.00} }, 
		{ {"NC", "CQ", "NC"}, { 70.0, 129.10} }, 
		{ {"O2", "C" , "O2"}, { 80.0, 126.00} }, 
		{ {"O2", "P" , "OH"}, { 45.0, 108.23} }, 
		{ {"OH", "P" , "OS"}, { 45.0, 102.60} }, 
		{ {"P" , "OS", "P" }, {100.0, 120.50} }, 
		{ {"C" , "CA", "HA"}, { 35.0, 120.00} },
		{ {"C" , "CB", "NB"}, { 70.0, 130.00} },
		{ {"C" , "CM", "CT"}, { 70.0, 119.70} },
		{ {"C" , "CM", "HA"}, { 35.0, 119.70} },
		{ {"C" , "CT", "H1"}, { 50.0, 109.50} },
		{ {"C" , "CT", "HP"}, { 50.0, 109.50} },
		{ {"C" , "CT", "N3"}, { 80.0, 111.20} },
		{ {"C" , "N" , "H" }, { 30.0, 120.00} },
		{ {"C" , "N*", "CT"}, { 70.0, 117.60} },
		{ {"C" , "NA", "C" }, { 70.0, 126.40} },
		{ {"C" , "NA", "H" }, { 30.0, 116.80} },
		{ {"C" , "OH", "HO"}, { 35.0, 113.00} },
		{ {"C*", "CB", "CN"}, { 63.0, 108.80} },
		{ {"C*", "CT", "HC"}, { 50.0, 109.50} },
		{ {"C*", "CW", "NA"}, { 70.0, 108.70} },
		{ {"CA", "C" , "OH"}, { 70.0, 120.00} },
		{ {"CA", "CA", "CB"}, { 63.0, 120.00} },
		{ {"CA", "CA", "CT"}, { 70.0, 120.00} },
		{ {"CA", "CA", "HA"}, { 35.0, 120.00} },
		{ {"CA", "CB", "CN"}, { 63.0, 116.20} },
		{ {"CA", "CM", "CM"}, { 63.0, 117.00} },
		{ {"CA", "CM", "HA"}, { 35.0, 123.30} },
		{ {"CA", "CN", "NA"}, { 70.0, 132.80} },
		{ {"CA", "CT", "HC"}, { 50.0, 109.50} },
		{ {"CA", "N2", "H" }, { 35.0, 120.00} },
		{ {"CA", "NC", "CB"}, { 70.0, 112.20} },
		{ {"CB", "C" , "NA"}, { 70.0, 111.30} },
		{ {"CB", "C*", "CT"}, { 70.0, 128.60} },
		{ {"CB", "CA", "H4"}, { 35.0, 120.00} },
		{ {"CB", "CA", "N2"}, { 70.0, 123.50} },
		{ {"CB", "CB", "N*"}, { 70.0, 106.20} },
		{ {"CB", "CB", "NC"}, { 70.0, 127.70} },
		{ {"CB", "N*", "CK"}, { 70.0, 105.40} },
		{ {"CB", "N*", "H" }, { 30.0, 125.80} },
		{ {"CB", "NC", "CQ"}, { 70.0, 111.00} },
		{ {"CC", "CT", "HC"}, { 50.0, 109.50} },
		{ {"CC", "CV", "NB"}, { 70.0, 120.00} },
		{ {"CC", "CW", "NA"}, { 70.0, 120.00} },
		{ {"CC", "NA", "H" }, { 30.0, 120.00} },
		{ {"CK", "N*", "CT"}, { 70.0, 128.80} },
		{ {"CM", "C" , "NA"}, { 70.0, 114.10} },
		{ {"CM", "CA", "N2"}, { 70.0, 120.10} },
		{ {"CM", "CM", "CT"}, { 70.0, 119.70} },
		{ {"CM", "CM", "HA"}, { 35.0, 119.70} },
		{ {"CM", "CT", "HC"}, { 50.0, 109.50} },
		{ {"CM", "N*", "H" }, { 30.0, 121.20} },
		{ {"CN", "NA", "CW"}, { 70.0, 111.60} },
		{ {"CR", "NA", "CW"}, { 70.0, 120.00} },
		{ {"CR", "NB", "CV"}, { 70.0, 117.00} },
		{ {"CT", "C" , "O" }, { 80.0, 120.40} },
		{ {"CT", "C*", "CW"}, { 70.0, 125.00} },
		{ {"CT", "CC", "CW"}, { 70.0, 120.00} },
		{ {"CT", "CC", "NB"}, { 70.0, 120.00} },
		{ {"CT", "CT", "H1"}, { 50.0, 109.50} },
		{ {"CT", "CT", "HC"}, { 50.0, 109.50} },
		{ {"CT", "CT", "N" }, { 80.0, 109.70} },
		{ {"CT", "CT", "N2"}, { 80.0, 111.20} },
		{ {"CT", "CT", "OH"}, { 50.0, 109.50} },
		{ {"CT", "CT", "S" }, { 50.0, 114.70} },
		{ {"CT", "N" , "CT"}, { 50.0, 118.00} },
		{ {"CT", "N2", "H" }, { 35.0, 118.40} },
		{ {"CT", "OH", "HO"}, { 55.0, 108.50} },
		{ {"CT", "OS", "P" }, {100.0, 120.50} },
		{ {"CT", "S" , "S" }, { 68.0, 103.70} },
		{ {"CV", "CC", "NA"}, { 70.0, 120.00} },
		{ {"CW", "CC", "NB"}, { 70.0, 120.00} },
		{ {"F" , "CT", "F" }, { 77.0, 109.10} },
		{ {"H" , "N" , "H" }, { 35.0, 120.00} },
		{ {"H" , "N3", "H" }, { 35.0, 109.50} },
		{ {"H1", "CT", "N" }, { 50.0, 109.50} },
		{ {"H1", "CT", "N2"}, { 50.0, 109.50} },
		{ {"H1", "CT", "OS"}, { 50.0, 109.50} },
		{ {"H1", "CT", "SH"}, { 50.0, 109.50} },
		{ {"H2", "CT", "N*"}, { 50.0, 109.50} },
		{ {"H4", "CM", "N*"}, { 35.0, 119.10} },
		{ {"H4", "CW", "NA"}, { 35.0, 120.00} },
		{ {"H5", "CK", "NB"}, { 35.0, 123.05} },
		{ {"H5", "CR", "NA"}, { 35.0, 120.00} },
		{ {"HC", "CT", "HC"}, { 35.0, 109.50} },
		{ {"HP", "CT", "HP"}, { 35.0, 109.50} },
		{ {"HS", "SH", "HS"}, { 35.0,  92.07} },
		{ {"N" , "C" , "O" }, { 80.0, 122.90} },
		{ {"N*", "C" , "NC"}, { 70.0, 118.60} },
		{ {"N*", "CB", "NC"}, { 70.0, 126.20} },
		{ {"N*", "CT", "OS"}, { 50.0, 109.50} },
		{ {"N2", "CA", "NA"}, { 70.0, 116.00} },
		{ {"NA", "C" , "O" }, { 80.0, 120.60} },
		{ {"NA", "CR", "NA"}, { 70.0, 120.00} },
		{ {"NC", "C" , "O" }, { 80.0, 122.50} },
		{ {"O" , "C" , "O" }, { 80.0, 126.00} },
		{ {"O2", "P" , "O2"}, {140.0, 119.90} },
		{ {"O2", "P" , "OS"}, {100.0, 108.23} },
		{ {"OS", "P" , "OS"}, { 45.0, 102.60} },
		{ {"X" , "X" , "X" }, {  0.0,   0.00} }
	};

	const static std::map<std::pair<String, String>, std::tuple<double, double, int, int>> s_defaultTorsion23Parameters = {
		{ {"C" , "CA"}, {14.50, 180.0, 2, 4} }, 
		{ {"C" , "CM"}, { 8.70, 180.0, 2, 4} }, 
		{ {"C" , "N" }, {10.00, 180.0, 2, 4} }, 
		{ {"C" , "NA"}, { 5.40, 180.0, 2, 4} }, 
		{ {"C" , "OH"}, { 1.80, 180.0, 2, 2} }, 
		{ {"C*", "CT"}, { 0.00,   0.0, 2, 6} }, 
		{ {"CA", "CA"}, {14.50, 180.0, 2, 4} }, 
		{ {"CA", "CM"}, {10.20, 180.0, 2, 4} }, 
		{ {"CA", "CT"}, { 0.00,   0.0, 2, 6} }, 
		{ {"CA", "NA"}, { 6.00, 180.0, 2, 2} }, 
		{ {"CB", "CB"}, {21.80, 180.0, 2, 4} }, 
		{ {"CB", "N*"}, { 6.60, 180.0, 2, 4} }, 
		{ {"CB", "NC"}, { 8.30, 180.0, 2, 2} }, 
		{ {"CC", "CV"}, {20.60, 180.0, 2, 4} }, 
		{ {"CC", "NA"}, { 5.60, 180.0, 2, 4} }, 
		{ {"CK", "N*"}, { 6.80, 180.0, 2, 4} }, 
		{ {"CM", "CM"}, {26.60, 180.0, 2, 4} }, 
		{ {"CM", "N*"}, { 7.40, 180.0, 2, 4} }, 
		{ {"CQ", "NC"}, {13.60, 180.0, 2, 2} }, 
		{ {"CR", "NB"}, {10.00, 180.0, 2, 2} }, 
		{ {"CT", "N" }, { 0.00,   0.0, 2, 6} }, 
		{ {"CT", "N2"}, { 0.00,   0.0, 3, 6} }, 
		{ {"CT", "OH"}, { 0.50,   0.0, 3, 3} }, 
		{ {"CT", "S" }, { 1.00,   0.0, 3, 3} }, 
		{ {"CV", "NB"}, { 4.80, 180.0, 2, 2} }, 
		{ {"OH", "P" }, { 0.75,   0.0, 3, 3} }, 
		{ {"C" , "CB"}, {12.00, 180.0, 2, 4} },
		{ {"C" , "CT"}, { 0.00,   0.0, 2, 4} },
		{ {"C" , "N*"}, { 5.80, 180.0, 2, 4} },
		{ {"C" , "NC"}, { 8.00, 180.0, 2, 2} },
		{ {"C*", "CB"}, { 6.70, 180.0, 2, 4} },
		{ {"C*", "CW"}, {26.10, 180.0, 2, 4} },
		{ {"CA", "CB"}, {14.00, 180.0, 2, 4} },
		{ {"CA", "CN"}, {14.50, 180.0, 2, 4} },
		{ {"CA", "N2"}, { 9.60, 180.0, 2, 4} },
		{ {"CA", "NC"}, { 9.60, 180.0, 2, 2} },
		{ {"CB", "CN"}, {12.00, 180.0, 2, 4} },
		{ {"CB", "NB"}, { 5.10, 180.0, 2, 2} },
		{ {"CC", "CT"}, { 0.00,   0.0, 2, 6} },
		{ {"CC", "CW"}, {21.50, 180.0, 2, 4} },
		{ {"CC", "NB"}, { 4.80, 180.0, 2, 2} },
		{ {"CK", "NB"}, {20.00, 180.0, 2, 2} },
		{ {"CM", "CT"}, { 0.00,   0.0, 3, 6} },
		{ {"CN", "NA"}, { 6.10, 180.0, 2, 4} },
		{ {"CR", "NA"}, { 9.30, 180.0, 2, 4} },
		{ {"CT", "CT"}, { 1.40,   0.0, 3, 9} },
		{ {"CT", "N*"}, { 0.00,   0.0, 2, 6} },
		{ {"CT", "N3"}, { 1.40,   0.0, 3, 9} },
		{ {"CT", "OS"}, { 1.15,   0.0, 3, 3} },
		{ {"CT", "SH"}, { 0.75,   0.0, 3, 3} },
		{ {"CW", "NA"}, { 6.00, 180.0, 2, 4} },
		{ {"OS", "P" }, { 4.80, 180.0, 2, 2} }
	};

	const static std::map<std::tuple<String, String, String, String>, std::vector<std::tuple<double, double, int, int>>> s_defaultTorsion1234Parameters = {
		{ {"C" , "N" , "CT", "C" }, { { 0.00,   0.0, -4, 1},   {0.00, 180.0, -3, 1},   {0.20, 180.0, -2, 1},   {0.00, 180.0,  1, 1} } },
		{ {"CT", "CT", "C" , "N" }, { { 0.10,   0.0, -4, 1},   {0.00,   0.0, -3, 1},   {0.07,   0.0, -2, 1},   {0.00, 180.0,  1, 1} } },
		{ {"CT", "CT", "N" , "C" }, { { 0.50, 180.0, -4, 1},   {0.15, 180.0, -3, 1},   {0.00, 180.0, -2, 1},   {0.53,   0.0,  1, 1} } },
		{ {"CT", "CT", "OS", "CT"}, { {0.383,   0.0, -3, 1},   {0.10, 180.0,  2, 1} } },
		{ {"CT", "S" , "S" , "CT"}, { {0.60,   0.0,  3,  1},   {3.50,   0.0, -2, 1} } },
		{ {"H" , "N" , "C" , "O" }, { {2.50, 180.0, -2,  1},   {2.00,   0.0,  1, 1} } },
		{ {"N" , "CT", "C" , "N" }, { {0.40, 180.0, -4,  1},   {0.00,   0.0, -3, 1},   {1.35, 180.0, -2, 1},   {0.75, 180.0,  1, 1} } },
		{ {"OH", "CT", "CT", "OH"}, { {0.144,   0.0, -3, 1},   {1.00,   0.0,  2, 1} } },
		{ {"OH", "P" , "OS", "CT"}, { {0.25,   0.0, -3,  1},   {1.20,   0.0,  2, 1} } },
		{ {"OS", "CT", "CT", "OH"}, { {0.144,   0.0, -3, 1},   {1.00,   0.0,  2, 1} } },
		{ {"OS", "CT", "CT", "OS"}, { {0.144,   0.0, -3, 1},   {1.00,   0.0,  3, 1} } },
		{ {"OS", "CT", "N*", "CK"}, { {0.50, 180.0, -2,  1},   {2.50,   0.0,  1, 1} } },
		{ {"OS", "CT", "N*", "CM"}, { {0.50, 180.0, -2,  1} } },
		{ {"OS", "P" , "OS", "CT"}, { {0.25,   0.0, -3,  1},   {1.20,   0.0, 2.0, 1} } },
		{ {"S" , "CT", "N*", "CM"}, { {2.50,   0.0,  1,  1} } }
	};

	const static std::map<std::pair<String, String>, double> s_defaultOutOfPlane34Parameters = {
		{ {"C" , "O" }, 10.5 },  
		{ {"CA", "HA"},  1.1 },  
		{ {"CM", "HA"},  1.1 },  
		{ {"CV", "H4"},  1.1 },  
		{ {"N2", "H" },  1.0 }, 
		{ {"CA", "H4"},  1.1 },
		{ {"CK", "H5"},  1.1 },
		{ {"CQ", "H5"},  1.1 },
		{ {"CW", "H4"},  1.1 },
		{ {"NA", "H" },  1.0 },
		{ {"CA", "H5"},  1.1 },
		{ {"CM", "H4"},  1.1 },
		{ {"CR", "H5"},  1.1 },
		{ {"N" , "H" },  1.0 }
	};

	const static std::map<std::tuple<String, String, String>, double> s_defaultOutOfPlane234Parameters = {
		{ { "CT", "N" , "CT" },  1.0 },
		{ { "N2", "CA", "N2" }, 10.5 },
		{ { "O2", "C" , "O2" }, 10.5 }
	};

	const static std::map<std::tuple<String, String, String, String>, double> s_defaultOutOfPlane1234Parameters = {
		{ {"CA", "CA", "C" , "OH"}, 1.1 }, 
		{ {"CB", "NC", "CA", "N2"}, 1.1 }, 
		{ {"CM", "C" , "CM", "CT"}, 1.1 }, 
		{ {"CT", "CM", "CM", "C" }, 1.1 }, 
		{ {"NC", "CM", "CA", "N2"}, 1.1 }, 
		{ {"NA", "CW", "CC", "CT"}, 1.1 }, 
		{ {"CA", "CA", "CA", "CT"}, 1.1 },
		{ {"CK", "CB", "N*", "CT"}, 1.0 },
		{ {"CM", "C" , "N*", "CT"}, 1.0 },
		{ {"CW", "CB", "C*", "CT"}, 1.1 },
		{ {"NA", "CV", "CC", "CT"}, 1.1 },
		{ {"NA", "NC", "CA", "N2"}, 1.1 },
		{ {"NB", "CW", "CC", "CT"}, 1.1 }
	};

	ForceField::ForceField() {
		SetDefaultParameters();
	}

	ForceField::ForceField(const String &filePath) {
		std::ifstream file(filePath);

		if (file.is_open()) {
			String line;

			while (std::getline(file, line)) {
				if (line.find("#") != 0 && !line.empty()) {
					String part = line.substr(0, line.find("#"));

					std::vector<String> tokens = utils::Tokenize(part);

					ResolveTokens(tokens);
				}
			}

			file.close();
		}
		else {
			std::cout << "Could not open force field file " << filePath << std::endl;
			std::cout << "Ignoring force field file and setting default force field parameters" << std::endl;
			SetDefaultParameters();
			return;
		}
	}

	std::ostream& operator<<(std::ostream &stream, const ForceField &forceField) {
		stream << "Force Field: {" << std::endl;

		stream << "\tAtomic masses: [" << std::endl;

		for (const auto& elem : forceField.m_atomicMasses) {
			stream << "\t\t(" << elem.first << ", " << elem.second << ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tElements radii: [" << std::endl;

		for (const auto& elem : forceField.m_covalentRadii) {
			stream << "\t\t(" << elem.first << ", " << elem.second << ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tVan Der Waals parameters: [" << std::endl;

		for (const auto& elem : forceField.m_vanDerWaalsParameters) {
			stream << "\t\t" << elem.first << ": (" << elem.second.first << ", " << elem.second.second <<  ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tBond length parameters: [" << std::endl;

		for (const auto& elem : forceField.m_bondLengthParameters) {
			stream << "\t\t(" << elem.first.first << ", " << elem.first.second << "): (" << elem.second.first << ", " << elem.second.second << ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tBond angle parameters: [" << std::endl;

		for (const auto& elem : forceField.m_bondAngleParameters) {
			stream << "\t\t(" << std::get<0>(elem.first) << ", " << std::get<1>(elem.second) << ", " << std::get<2>(elem.first) << "): (";
			stream << elem.second.first << ", " << elem.second.second << ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tTwo-atom torsion parameters: [" << std::endl;

		for (const auto& elem : forceField.m_torsion23Parameters) {
			stream << "\t\t(" << elem.first.first << ", " << elem.first.second << "): (";
			stream << std::get<0>(elem.second) << ", " << std::get<1>(elem.second) << ", " << std::get<2>(elem.second) << ", " << std::get<3>(elem.second) << ")" << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tFour-atom torsion parameters: {" << std::endl;

		for (const auto& elem : forceField.m_torsion1234Parameters) {
			stream << "\t\t(" << std::get<0>(elem.first) << ", " << std::get<1>(elem.first) << ", " << std::get<2>(elem.first) << ", " << std::get<3>(elem.first) << "): [" << std::endl;
			
			for (const auto &tuple : elem.second) {
				stream << "\t\t\t(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ", " << std::get<2>(tuple) << ", " << std::get<3>(tuple) << ")" << std::endl;
			}

			stream << "\t\t]" << std::endl;
		}

		stream << "\t}" << std::endl;

		stream << "\tTwo-atom out-of-plane parameters: [" << std::endl;

		for (const auto& elem : forceField.m_outOfPlane34Parameters) {
			stream << "\t\t(" << elem.first.first << ", " << elem.first.second << "): " << elem.second << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tThree-atom out-of-plane parameters: [" << std::endl;

		for (const auto& elem : forceField.m_outOfPlane234Parameters) {
			stream << "\t\t(" << std::get<0>(elem.first) << ", " << std::get<1>(elem.first) << ", " << std::get<2>(elem.first) << "): " << elem.second << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "\tFour-atom out-of-plane parameters: [" << std::endl;

		for (const auto& elem : forceField.m_outOfPlane1234Parameters) {
			stream << "\t\t(" << std::get<0>(elem.first) << ", " << std::get<1>(elem.first) << ", " << std::get<2>(elem.first) << ", " << std::get<3>(elem.first) << "): ";
			stream << elem.second << std::endl;
		}

		stream << "\t]" << std::endl;

		stream << "}";

		return stream;
	}

	void ForceField::SetDefaultParameters() {
		m_atomicMasses = s_defaultAtomicMasses;
		m_covalentRadii = s_defaultCovalentRadii;
		m_vanDerWaalsParameters = s_defaultVanDerWaalsParameters;
		m_bondLengthParameters = s_defaultBondLengthParameters;
		m_bondAngleParameters = s_defaultBondAngleParameters;
		m_torsion23Parameters = s_defaultTorsion23Parameters;
		m_torsion1234Parameters = s_defaultTorsion1234Parameters;
		m_outOfPlane34Parameters = s_defaultOutOfPlane34Parameters;
		m_outOfPlane234Parameters = s_defaultOutOfPlane234Parameters;
		m_outOfPlane1234Parameters = s_defaultOutOfPlane1234Parameters;
	}

	void ForceField::ResolveTokens(const std::vector<String> &tokens) {
		if (tokens[0] == "mass" && tokens.size() >= 3) {
			m_atomicMasses.insert(std::make_pair(tokens[1], utils::ToDouble(tokens[2])));
		}

		if (tokens[0] == "radius" && tokens.size() >= 3) {
			m_covalentRadii.insert(std::make_pair(tokens[1], utils::ToDouble(tokens[2])));
		}

		if (tokens[0] == "vdw" && tokens.size() >= 3) {
			m_vanDerWaalsParameters.insert(std::make_pair(tokens[1], std::make_pair(utils::ToDouble(tokens[2]), utils::ToDouble(tokens[3]))));
		}

		if (tokens[0] == "bond" && tokens.size() >= 4) {
			m_bondLengthParameters.insert(std::make_pair(std::make_pair(tokens[1], tokens[2]), 
				std::make_pair(utils::ToDouble(tokens[3]), utils::ToDouble(tokens[4]))));
		}

		if (tokens[0] == "angle" && tokens.size() >= 6) {
			m_bondAngleParameters.insert(std::make_pair(std::make_tuple(tokens[1], tokens[2], tokens[3]), 
				std::make_pair(utils::ToDouble(tokens[4]), utils::ToDouble(tokens[5]))));
		}

		if (tokens[0] == "torsion23" && tokens.size() >= 6) {
			m_torsion23Parameters.insert(std::make_pair(std::make_pair(tokens[1], tokens[2]), 
				std::make_tuple(utils::ToDouble(tokens[3]), utils::ToDouble(tokens[4]), utils::NextInt(tokens[5]), utils::NextInt(tokens[6]))));
		}

		if (tokens[0] == "torsion1234" && (tokens.size() - 5) > 0 && (tokens.size() - 5) % 4 == 0) {
			std::pair<std::tuple<String, String, String, String>, std::vector<std::tuple<double, double, int, int>>> torsion1234Parameter;

			torsion1234Parameter.first = std::make_tuple(tokens[1], tokens[2], tokens[3], tokens[4]);

			int i = 5;
			
			while (i < tokens.size()) {
				torsion1234Parameter.second.push_back(std::make_tuple(
					utils::ToDouble(tokens[i]), utils::ToDouble(tokens[i + 1]), utils::NextInt(tokens[i + 2]), utils::NextInt(tokens[i + 3])));

				i += 4;
			}

			m_torsion1234Parameters.insert(torsion1234Parameter);
		}

		if (tokens[0] == "outOfPlane34" && tokens.size() >= 4) {
			m_outOfPlane34Parameters.insert(std::make_pair(std::make_pair(tokens[1], tokens[2]), utils::ToDouble(tokens[3])));
		}

		if (tokens[0] == "outOfPlane234" && tokens.size() >= 5) {
			m_outOfPlane234Parameters.insert(std::make_pair(std::make_tuple(tokens[1], tokens[2], tokens[3]), utils::ToDouble(tokens[4])));
		}

		if (tokens[0] == "outOfPlane1234" && tokens.size() >= 6) {
			m_outOfPlane1234Parameters.insert(std::make_pair(std::make_tuple(tokens[1], tokens[2], tokens[3], tokens[4]), utils::ToDouble(tokens[5])));
		}
	}

}