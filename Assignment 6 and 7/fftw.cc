#include "fftw.h"

// This is the function that computes the FFTW of the signal.

rarray<std::complex<double>,1> FFTW(rarray<std::complex<double>,1> signal) {

	// Create an array to store the results.
	rarray<std::complex<double>,1> f_t(signal.size());

	// Compute the FFtW.
	fftw_plan plan = fftw_plan_dft_1d(signal.size(), (fftw_complex*)signal.data(), (fftw_complex*)f_t.data(), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	// Return our computed Fourier transform.
	return f_t;
}
