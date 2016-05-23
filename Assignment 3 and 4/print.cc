// print.cc
// This is the module that prints the time step and the total number of ants on the table at that time step.
#include "print.h"
using namespace std; // Removes need of using prefix 'std'.

void Print(int t, float &totants) {
    cout << t << " " << totants << std::endl; // Output per step.
}
