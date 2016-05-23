#define BOOST_TEST_MODULE ants_test // Add our pseudo preprocessor guards here.
#define BOOST_TEST_DYN_LINK // Define here so no need to include flag when compiling.
#include <boost/test/included/unit_test.hpp> // Does not work without 'included' part.
#include "single_time_step.h"
#include "initial_speed.h"
#include "initial_amount.h"

// The test essentially mimics the ants driver code without looping through a series of time steps. It will test to see if the number of ants left after a time step is less than the amount before the time step, as is expected.

BOOST_AUTO_TEST_SUITE(Ants_test)
BOOST_AUTO_TEST_CASE(Ants_test) 
{
    const float Array_Size_X = 356;
    const float Array_Size_Y = 356;
    const int T_Start = 0; // Start time of simulation.
    const int T_End = 40; // End time of simulation.
    const int total_ants = 1010;
    const int Initial_I = 0; // Initial value of i.
    const int Initial_J = 0; // Initial value of j.
    float new_totants;

    rarray<float,2> number_of_ants(Array_Size_X,Array_Size_Y);
    rarray<float,2> velocity_of_ants(Array_Size_X,Array_Size_Y);
    rarray<float,2> new_number_of_ants(Array_Size_X,Array_Size_Y); // We must initialize this array here as it is called for within the function.

    Initial_Speed(velocity_of_ants); // Sets the velocity.

    Initial_Amount(number_of_ants, total_ants); // Sets the initial number.

    float totants;
    totants = 0.0;

    Single_Time_Step(number_of_ants, velocity_of_ants, new_number_of_ants, totants); // This is the individual time step process that we are testing without the for loop.
    // Now we will sum a new total number of ants to compare to our original.
    for (int i=Initial_I;i<Array_Size_X;i++) {
        for (int j=Initial_J;j<Array_Size_Y;j++) {
            new_totants += number_of_ants[i][j];
        }
    }
    BOOST_CHECK(new_totants < total_ants); // This is the piece de resistance, the actual check.
}
BOOST_AUTO_TEST_SUITE_END ()
// Compiles with no errors detected!



