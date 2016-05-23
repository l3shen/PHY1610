// simulation.cc
// This module describes what occurs during each time step, up to a final time.
#include "single_time_step.h"

using namespace std; // Removes need to prefix with std.

void Single_Time_Step(rarray<float,2> &number_of_ants, rarray<float,2> &velocity_of_ants, rarray<float,2> new_number_of_ants, float &totants) { // t is called for print function to work properly.
    // These are our constants used within our module.
    const int Initial_I = 0; // Initial value of i.
    const int Initial_J = 0; // Initial value of j.
    const float Array_Size_X = 356; // Initial X value for array size.
    const float Array_Size_Y = 356; // Initial Y value for array size.

    // rarray<float,2> new_number_of_ants(Array_Size_X,Array_Size_Y); // We must initialize this array here as it is called for within the function.

    for (int i=Initial_I;i<Array_Size_X;i++) {
        for (int j=Initial_J;j<Array_Size_Y;j++) {
            totants += number_of_ants[i][j];
        }
    }
    
    // Print(t, totants); // Print before doing the math, such that the first value is not the sum of the initial plus the first time step.

    new_number_of_ants.fill(0.0); // Replace loop with rarray fill method. 

    for (int i=Initial_I;i<Array_Size_X;i++) { // What happens per step.
            // Some ants fall, some remain on the table.
        for (int j=Initial_J;j<Array_Size_Y;j++) {
            int di = 1.9*sin(velocity_of_ants[i][j]);
            int dj = 1.9*cos(velocity_of_ants[i][j]); 
            int i2 = i + di;
            int j2 = j + dj;
                // Some ants do not walk on the table.
            new_number_of_ants[i][j]+=0.8*number_of_ants[i][j];
                // The remainder will walk and some will fall off.
            if (i2>=0 and i2<Array_Size_X and j2>=0 and j2<Array_Size_Y) { // Bug amended here. Previously strange restrictions.
                new_number_of_ants[i2][j2]+=0.2*number_of_ants[i][j];
            }
        }
    }
    number_of_ants = new_number_of_ants.copy(); // Replace loop with rarray copy method.
}
