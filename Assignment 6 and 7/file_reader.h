// file_reader.h
#ifndef FILE_READER_H
#define FILE_READER_H
#include <fstream>
#include <iostream>
#include <rarray>
#include <rarrayio>
#include <string>
#include <complex>
#include <cmath>
void Reading(const std::string filename, rarray<double,1> &times, rarray<std::complex<double>,1> &signal);
#endif
