#include "power_spectrum_cc.h"

// This function will calculate the power spectrum.
rarray<double,1> Power_Spectrum_Calc(rarray<std::complex<double>,1> FT_Signal) {

	// Parameters for the for loop.
	int start = 0;
	int end = (FT_Signal.size()-1);

	// We will store our results and output in an rarray.
	rarray<double,1> result(FT_Signal.size());

	// The for loop will compute the power spectrum.
	for (int n = start; n <= end; n++) {
		result[n] = real(FT_Signal[n]*std::conj(FT_Signal[n]));
	}

	return result;
}

// This function calculates the CC as defined within the assignment.
double Correlation_Coefficient(rarray<double,1> Power_Spectrum, rarray<double,1> Power_Spectrum_Det) {

	// Define our variable to be returned.
	double CC;

	// Now we calculate our CC.
	double Correlation_1 = cblas_ddot(Power_Spectrum.size(), Power_Spectrum.data(), 1, Power_Spectrum_Det.data(), 1);
	// Calculating the two portions of the denominator.
	double Correlation_2 = cblas_ddot(Power_Spectrum.size(), Power_Spectrum.data(), 1, Power_Spectrum.data(), 1);
	double Correlation_3 = cblas_ddot(Power_Spectrum.size(), Power_Spectrum_Det.data(), 1, Power_Spectrum_Det.data(), 1);
	// Calculating the final correlation coefficient.
	CC = (Correlation_1)/(std::sqrt(Correlation_2*Correlation_3));

	return CC;
}
