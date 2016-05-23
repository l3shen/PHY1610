// fftw.h
#ifndef FFTW_H
#define FFTW_H
#include <cmath>
#include <fftw3.h>
#include <rarray>
#include <complex>
#include <iostream>
rarray<std::complex<double>,1> FFTW(rarray<std::complex<double>,1> signal);
#endif
