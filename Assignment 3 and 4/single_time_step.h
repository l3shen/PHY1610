// single_time_step.h
// This is the time step simulation portion of the overall code.
// This part describes what occurs per time step to an end time. 

#ifndef SIMULATION_H
#define SIMULATION_H
#include <rarray>
#include <cmath>
#include <fstream> // As we print out a file (just for funsies).
#include "print.h"
#include "ticktock.h"
void Single_Time_Step(rarray<float,2> &number_of_ants, rarray<float,2> &velocity_of_ants, rarray<float,2> new_number_of_ants, float &totants);
#endif
