// Time-step module for Monte Carlo.

#include "walkring_timestep.h"
#include <random>

using namespace std;

// Initialization for RNG stuff.
default_random_engine engine; 
uniform_real_distribution<double> uniform(0.0,1.0); // Use 'uniform' to ensure equally probably random numbers.

// The probabilities are outlined in the assignment and have been applied here. They are calculated in the main driver.

void perform_time_step(rarray<int,1>& pos, int N, double p)
{
	// Determine the total number of walkers from the matrix.
  	int number_of_walkers = pos.size() - 1;

	// Iterate over all walkers, which will generate a new random number to compare with the probability over each walker. 
  	for (int i = 0; i <= number_of_walkers; i++) {

		// Initialize our random double for determing probabilities. 
    		double r = uniform(engine);

		// Determine if the walker moves to the left.
    		if (r > 0 and r <= p) {
      			// Moving left.
      			pos[i] -= 1;
      			// If the new position is negative, set the new position at the other 'edge' of the ring.
      			if (pos[i] < 0) {
				pos[i] = N - 1;
      			}
    		}

		// Determine if the walker moves to the right.
    		else if (r > p and r <= (2*p)) {
      			// Moving right.
      			pos[i] += 1;
      			// boundary periodicity
      			if (pos[i] > N - 1) {
				pos[i] = 0;
      			}
    		}	
	
		// There is no need to determine if it doesn't move as if the randomly generated number does not fit the given criteria it will not move.

  	}
}


