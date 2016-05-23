// This is our driver code.
#include <iostream> 
#include <iomanip> // Added to be able to set precision of output values, for cleanliness.
#include <fstream> // Needed to print to file.
#include <rarray> // Used to store data from the f_all_min function in the form of an array, to export more than one value out of the function.
#include "f_all_min.h" // Header of the minima finding function file.

using namespace std; // Removes a need to write 'std' over and over.

int main() {
	rarray<double,2>Results(2,2); // Initialize our array to carry our minima and energy values.

	// We will initialize our variables outside of the for loop to make processing a little less intensive.
	double x1 = 0; // Our first local minima position.
	double x2 = 0; // Our second local minima position.
	double en1 = 0; // Our first local minima energy.
	double en2 = 0; // Our second loca minima energy.

	ofstream fout("global_minima.txt"); // Name of the output file.

	// Introducing our output first.
	cout << "First, here we have the given and mass and the x-position of the global minimum." << endl;
	cout << "The following information may be found in the text file 'global_minima.txt'." << endl;

	// Printing our column headers.
	cout << "Mass (kg)" << "\t" << "Minimum (x)" << endl;
	fout << "Mass (kg)" << "\t" << "Minimum (x)" << endl;

	// Begin looping over mass values 0 to 0.5 kg with 25 iterations.
	// Determining our range:
	double Minimum = 0.0;
	double Maximum = 0.5;
	double Points = 25;
	double Increment_of_Mass = (Maximum - Minimum)/(Points - 1);	

	for (double m = 0.00; m <= .50; m += Increment_of_Mass) {
 		Results = f_all_min(m); // Compute our results in to an array.
	
		// Set our variables to the energy and x positions from the array.
		x1 = Results[0][0];
		x2 = Results[0][1];
		en1 = Results[1][0];
		en2 = Results[1][1];

		// We will now compare the two to find the global minima. The program will only print the lowest energy x position, i.e. the global minima.
		// This is accomplished with the following if/else statement.
		if (en1 < en2) {
			cout << setprecision(3) << m << fixed << "\t" << setprecision(5) << x1 << fixed << endl;
			fout << setprecision(3) << m << fixed << "\t" << setprecision(5) << x1 << fixed << endl;
		}
		else {
			cout << setprecision(3) << m << fixed << "\t" << setprecision(5) << x2 << fixed << endl;
			fout << setprecision(3) << m << fixed << "\t" << setprecision(5) << x2 << fixed << endl;
		}
	}

	// Close our file output and return a value to the shell.
	fout.close(); 
	
	// This portion of the file will use the root finder to find a minima in the value of the difference of potentials, effectively giving us the maximum load.
	// We set up our range (cannot use 0 or 0.5 as they are infinite).
	double x_lo = 0.001;
	double x_hi = 0.499;
	// Initialize and run our solver.
	gsl_root_fsolver *solver;
	gsl_function fwrapper;
	solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
	fwrapper.function = energy_difference;
	fwrapper.params = 0; // We are not using any parameters here. 
	gsl_root_fsolver_set(solver, &fwrapper, x_lo, x_hi);
	int status = 1;
	double x_rt;
	for (int iter = 0; status and iter < 100; ++iter) {
		gsl_root_fsolver_iterate(solver);
		x_rt = gsl_root_fsolver_root(solver);
		double x_lo = gsl_root_fsolver_x_lower(solver);
		double x_hi = gsl_root_fsolver_x_upper(solver);
		status = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
	}
	gsl_root_fsolver_free(solver);	
	cout << "The maximum load is at: " << x_rt << " m" << endl;

	return status;

}
