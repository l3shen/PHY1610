// This is the driver script. It links all of the modules together.
#include <cmath>
#include <iostream>
#include <rarray> // Including the rarray library.
// Here we link the header files. Not that the print header file is not included, as that module is called by simulation.cc and not by the main driver.
#include "initial_speed.h"
#include "initial_amount.h"
#include "single_time_step.h" // Renamed from  simulation.
#include "ticktock.h" // Adding ticktock library.

// The file was split in to four modules: two for initialization, one for the actual simulation, and print, which is called in the simulation.

// This is our main driver file that puts everything together.
// Note that each function pass a value by reference; this is important!
int main()
{
    TickTock stopwatch; // Initialize our ticktock library.
    stopwatch.tick(); // Begin timing our simulation.
    // The only constant we need in our main driver loop are the array size. 
    // I adjusted the constant to two variables, such that if somebody wanted to create a rectangular array it is possible to do so.
    const float Array_Size_X = 356;
    const float Array_Size_Y = 356;

    // Here are the relevant constants for the time steps. Can adjust length of time period.
    const int T_Start = 0; // Start time of simulation.
    const int T_End = 40; // End time of simulation.

    const int total_ants = 1010;

    // Finally, here is where the simulation occurs.
    // First, we define our arrays. They must be defined in order for the functions below to work (i.e. know what to use when called).
    // Notice how the new_number_of_ants array was not needed; this is because the functions here do not call for it, although it is called in the simulation function and thus initialized there.
    rarray<float,2> number_of_ants(Array_Size_X,Array_Size_Y);
    rarray<float,2> velocity_of_ants(Array_Size_X,Array_Size_Y);
    rarray<float,2> new_number_of_ants(Array_Size_X,Array_Size_Y); // We must initialize this array here as it is called for within the function.

    // Initialization of ants.
    // First we define the speed of the ants.
    stopwatch.tick(); // Added timer for initialization (speed).
    Initial_Speed(velocity_of_ants); // Sets the velocity.
    stopwatch.tock("Time for initialization of time:");

    stopwatch.tick(); // Added timer for initialization (number).
    Initial_Amount(number_of_ants, total_ants); // Sets the initial number.
    stopwatch.tock("Time for initialization of number:");

    // Simulation time step and print out.
    for (int t = T_Start; t < T_End; t++) { // Begin our overall for loop.
        float totants;
        totants = 0.0;
        Single_Time_Step(number_of_ants, velocity_of_ants, new_number_of_ants, totants); // This is the individual time step process. It is now modularized properly and does not contain the foor loop - it is a true single time step module now.
        Print(t, totants); // Since print is no longer called for in the old module, it may be placed in the main driver file.
    }
    stopwatch.tock("The total simulation time is:");
    // Perhaps delete rarrays afterwards to clear memory? I did not include that, but it may be a possible improvement.
    return 0; // Return a value to the shell.
}             


