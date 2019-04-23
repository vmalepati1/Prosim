#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

/* Conversion of acceleration from [kcal/(A*g) to [A/(ps^2)]. */
#define ACCELERATION_CONVERSION 418.400000

/* Threshold beyond covalent radii sum to determine bond cutoff. */
#define BOND_THRESHOLD 1.2
/* Displacement distance [Angstrom] for numerical gradient. */
#define NUMERICAL_DISPLACEMENT 1.0E-6

#define RADIANS_TO_DEGREES 180.0 / M_PI
#define DEGREES_TO_RADIANS M_PI / 180.0

/* Conversion of electrostatic energy from[ceu] to[kcal / mol]. */
#define CEU_TO_KCAL 332.06375
/* Conversion of kinetic energy from[amu*A ^ 2 / ps ^ 2] to[kcal / mol]. */
#define KINETIC_TO_KCAL 0.00239005736
/* Conversion from [kcal*A^3/mol] to [Pa] for pressure. */
#define KCAL_A_MOL_TO_PA 69476.95

/* Boltzmann constant [kcal/(mol*K)]. */
#define BOLTZMANN_CONSTANT 0.001987204

/* Gas constant in units of [amu*A^2/(ps^2*K)]. */
#define GAS_CONSTANT 0.83144598