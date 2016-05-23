// pot_energy.h

#ifndef POT_ENERGY_H
#define POT_ENERGY_H
#include <cmath>
// We will define our struct here so we do not have to do it everywhere else.
struct Params{
	double a, b, c, d, f, g, m;
};
double pot_energy(double x, void* param);
#endif
