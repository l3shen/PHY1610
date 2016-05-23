// power_spectrum_cc.h
#ifndef POWER_SPECTRUM_CALC_H
#define POWER_SPECTRUM_CALC_H
#include <cblas.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <rarray>
#include <rarrayio>
rarray<double,1> Power_Spectrum_Calc(rarray<std::complex<double>,1> FT_Signal);
double Correlation_Coefficient(rarray<double,1> Power_Spectrum, rarray<double,1> Power_Spectrum_Det);
#endif
