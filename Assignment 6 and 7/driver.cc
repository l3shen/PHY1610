// This is the driver file that puts it all together.
// Kamil Krawczyk - PHY1610

#include "file_reader.h"
#include "fftw.h"
#include "filename.h"
#include "power_spectrum_cc.h"
#include "highest_value.h"

int main() {

	// Print our header to console.
	std::cout << "Welcome to the gravitational waves correlation coefficient calculator!" << '\n' << \
	"The detected signal file number and the correlation coefficient to the predicted signal:" << '\n' << std::endl;

	// Initialize our rarrays. This is for the predicted signal.
	rarray<double,1> times;
	rarray<std::complex<double>,1> signal;

	// Reading data from the .rat files for the predicted signal.
	Reading("GWprediction.rat", times, signal);
	
	// Creating an array to house the FT data for the predicted signal. Calculate the power spectrum of the predicted signal.
	rarray<std::complex<double>,1> FT_Signal = FFTW(signal);
	rarray<double,1> Power_Spectrum(signal.size());
	Power_Spectrum = Power_Spectrum_Calc(FT_Signal);

	// Again, initialize rarrays. For the detected signal.
	rarray<double,1> times_det;
	rarray<std::complex<double>,1> signal_det;

	// Define an array to hold the value of the coefficients upon completion of the for loop.
	rarray<double,1> Coefficients(32); // Defined with a size that corresponds to the number of data points.

	// Initialize parameters for the maximimum finding part.
	int N_Max = 5; // The number of sorted values we need.
	rarray<int,1> File_Index(N_Max); // A placeholder array that will hold the index of the file name so we can keep track; passed by reference in to our function.


	// Here is a loop that iterates through the inidividual data files, computes their power spectrum, then does the correlation computation and dumps it in to an rarray.
	for (int z = 1; z <= 32; z++) {

		// Define the filename of the detected signal.
		const std::string file = Filename(z);

		// Reading data from the .rat files of the detected signal. 
		Reading(file, times_det, signal_det);

		// Creating an array to house the FT data for the detected signal.
		rarray<std::complex<double>,1> FT_Signal_Det = FFTW(signal_det);
		rarray<double,1> Power_Spectrum_Det(signal_det.size());	
		Power_Spectrum_Det = Power_Spectrum_Calc(FT_Signal_Det);

		double CC = Correlation_Coefficient(Power_Spectrum, Power_Spectrum_Det);
	
		std::cout << "File number: " << z << ". The correlation coefficient: " << CC << std::endl;	

		Coefficients[z] = CC; // Insert the coefficient values in to our previously defined array for later.

		}

	// Compute the highest values within the array of correlation coefficients.
	rarray<double,1> HighestValues = Highest_Value(Coefficients, File_Index, N_Max);

	// Print another header to the console.
	std::cout << '\n' << "The top five correlations, their file index, and their correlation coefficients:" << std::endl;
	
	// Print the results and the corresponding file indices.
	Print(File_Index, HighestValues, N_Max);

	return 0; // Return a value to the shell.
}

