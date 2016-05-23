#include "f_all_min.h"
#include "pot_energy.h"

// Since our function contains either one or two minima, we will initialize two solvers to find both minima. It it set up in such a way that if there is only one minimum,
// it will simply return two of the same, which then can be compared using a for loop in the driver code. 

// The function itself returns an rarray because that was the easiest means of outputting more than one value from this function; initially, I had used
// a double type for the function but alas the caveat was that I could only return one value. Hence the use of rarray.
using namespace std;

rarray<double,2> f_all_min(double mass) {
	// We will define our search range and initial guess for the first minimum here. We will also define our tolerance.
	double tolerance = 0.0001;
	double x_lo_1 = 0.01;
	double x_hi_1 = 0.4999;
	double m1 = 0.14;

	// We will initialize our array for storing our data here.
	rarray<double, 2> Values(2,2);

	// We will initialize our solver here.
	// Our parameters here can be changed at any time; I chose not to hard code them to make the function more "moddable".
	Params args = {1, 0.1, 100, 0.5, 2500, 9.8, mass};
	
	// We will initialize some pointers and select our solver.
	gsl_min_fminimizer* solver;
	gsl_function fwrapper;
	solver = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
	fwrapper.function = pot_energy;
	fwrapper.params = &args;

	// We will now set the requisite parameters for our solver - the solver type, the function, and our parameters.
	gsl_min_fminimizer_set(solver, &fwrapper, m1, x_lo_1, x_hi_1);
	int status1 = 1;
	double x_min_1;
	for (int iter=0; status1 and iter < 100; ++iter) {
		gsl_min_fminimizer_iterate(solver);
		x_min_1 = gsl_min_fminimizer_x_minimum(solver);
		double x_lo_1 = gsl_min_fminimizer_x_lower(solver);
		double x_hi_1 = gsl_min_fminimizer_x_upper(solver);
		status1 = gsl_min_test_interval(x_lo_1,x_hi_1,0,tolerance);
	}

	// We will now calculate the energy at this minima.
	double Energy1 = pot_energy(x_min_1, &args);

	// We will now place the results in our array.
	Values[0][0] = x_min_1;
	Values[1][0] = Energy1;

	// Simiarly, we will set up another search.
	double x_lo_2 = 0.4;
	double x_hi_2 = 0.4999;
	double m2 = 0.42;

	gsl_set_error_handler_off (); // We will put in our own error check.

	gsl_min_fminimizer_set(solver, &fwrapper, m2, x_lo_2, x_hi_2);
	int status = 1;
	double x_min_2;
	for (int iter=0; status and iter < 100; ++iter) {
		gsl_min_fminimizer_iterate(solver);
		x_min_2 = gsl_min_fminimizer_x_minimum(solver);
		double x_lo_2 = gsl_min_fminimizer_x_lower(solver);
		double x_hi_2 = gsl_min_fminimizer_x_upper(solver);
		status = gsl_min_test_interval(x_lo_2,x_hi_2,0,tolerance);
	}
	gsl_min_fminimizer_free(solver);

	double Energy2 = pot_energy(x_min_2, &args);

	// Here will check for the validity of the second minima.
	if (pot_energy(x_min_2*(1+tolerance),&args) > Energy2 and pot_energy(x_min_2*(1-tolerance),&args) > Energy2) {
		Values[0][1] = x_min_2;
		Values[1][1] = Energy2;
	}
	else {
		cout << "Error: There is only one minimum." << endl;
	}

	return Values;
	
}

// This is our function to compare our energies; when this is at '0', i.e. the root, it is at the maximum extension.
double energy_difference(double mass, void * param) {
	rarray<double, 2> Values_Comparison(2,2);
	
	Values_Comparison = f_all_min(mass);

	double First_Energy = Values_Comparison[1][0];
	double Second_Energy = Values_Comparison[1][1];

	double Energy_Difference = Second_Energy - First_Energy;

	return Energy_Difference;
}

		
