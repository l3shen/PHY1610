#include "pot_energy.h"

// Here is our function for the total potential energy. Both terms have been combined for simplicity; otherwise, I would just make a total of three functions.
double pot_energy(double x, void* param) {
	Params* p = (Params*)param;
	return (p->a*((p->b/x)+((pow(p->d,2))/(p->f*(pow((x-p->d),2)))-exp((-p->c*((pow((x-p->b),2))/(2*p->a))))))) + (-p->m*p->g*x);
}
