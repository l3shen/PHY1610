#include "file_reader.h"

// This function reads in the file, as per the assignment notes!
void Reading(const std::string filename, rarray<double,1> &times, rarray<std::complex<double>,1> &signal) {
	std::ifstream f(filename.c_str());
	// Read in the signal.
	f >> times;
	f >> signal;
}
