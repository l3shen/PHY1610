// initial_amount.cc
// This is the second module declaring the initial amount of ants.
#include "initial_amount.h"

void Initial_Amount(rarray<float,2> &number_of_ants, const int total_ants) {
    // Constants used. Again, included in the function not the header, although the opposite may be done.
    const float Array_Size_X = 356; // Initial X value for array size.
    const float Array_Size_Y = 356; // Initial Y value for array size.
    const int Initial_I = 0; // Initial value of i.
    const int Initial_J = 0; // Initial value of j.
    const int Start_N = 0; // The counting value for ants.
    const float Start_Z = 0; // Not sure what this is for for now.
    // const int total_ants = 1010; // Initial number of ants.
    // Defining our n and z values outside of the for loop.
    int n = Start_N;
    float z = Start_Z;
    // Function.
    for (int i=Initial_I;i<Array_Size_X;i++) {
        for (int j=Initial_J;j<Array_Size_Y;j++) {
            number_of_ants[i][j] = 0.0; // Method within rarray array_name.fill(value)
        }
    }
    while (n < total_ants) {
        for (int i=Initial_I;i<Array_Size_X;i++) {
            for (int j=Initial_J;j<Array_Size_Y;j++) {
                z += sin(0.3*(i+j));
                if (z>1 and n!=total_ants) {
                    number_of_ants[i][j] += 1;
                    n += 1;
                }
            }
        }
    }
}
