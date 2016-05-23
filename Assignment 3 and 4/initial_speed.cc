// initial_speed.cc
// This is our module for the initialization of the speed of the ants.
#include "initial_speed.h"

void Initial_Speed(rarray<float,2> &velocity_of_ants) {
     // Constants used in the function. These can be included in the header but I left them in their respective function codes instead.
     const float Array_Size_X = 356; // Initial X value for array size.
     const float Array_Size_Y = 356; // Initial Y value for array size.
     const int Initial_I = 0; // Initial value of i.
     const int Initial_J = 0; // Initial value of j.
     const float Speed_Modifier = 3560; // Constant to modify initial velocity.
     // The function itself.
     for (int i=Initial_I;i<Array_Size_X;i++) {
        for (int j=Initial_J;j<Array_Size_Y;j++) {
            velocity_of_ants[i][j] = 3.1415*(sin((2*3.1415*(i+j))/Speed_Modifier)+1);
            // The velocity of ants can be converted to a module in itself, but I chose to not overly complicate things for the time being.
        }
    }
}
