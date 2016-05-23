// f_all_min.h
#ifndef F_ALL_MIN_H
#define F_ALL_MIN_H
#include <iostream>
#include <cfloat>
#include <utility>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include "pot_energy.h"
#include <rarray>
rarray<double,2> f_all_min(double m);
double energy_difference(double mass, void * param);
#endif
